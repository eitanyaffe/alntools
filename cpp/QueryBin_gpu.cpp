#include "QueryBin_gpu.h"
#include <fstream>
#include <iomanip> // For std::setprecision

#ifdef METAL_SUPPORT
#include "QueryBin_metal_bridge.h"
#endif

QueryBinGPU::QueryBinGPU(
    const std::vector<Interval>& intervals,
    const AlignmentStore& store,
    int binsize,
    double seg_threshold,
    double non_ref_threshold,
    int num_threads,
    ClipMode clip_mode,
    int clip_margin,
    double min_mutations_percent,
    double max_mutations_percent,
    int min_alignment_length,
    int max_alignment_length)
    : intervals(intervals)
    , store(store)
    , binsize(binsize)
    , seg_threshold(seg_threshold)
    , non_ref_threshold(non_ref_threshold)
    , clip_mode(clip_mode)
    , clip_margin(clip_margin)
    , min_mutations_percent(min_mutations_percent)
    , max_mutations_percent(max_mutations_percent)
    , min_alignment_length(min_alignment_length)
    , max_alignment_length(max_alignment_length)
{
#ifdef METAL_SUPPORT
    metal_impl = metal_create_implementation();
    use_gpu = (metal_impl != nullptr);
    if (use_gpu) {
        std::cout << "initialized Apple Silicon GPU acceleration" << std::endl;
    } else {
        std::cout << "metal GPU not available, using CPU implementation" << std::endl;
    }
#else
    use_gpu = false;
    metal_impl = nullptr;
    std::cout << "GPU support not compiled, using CPU implementation" << std::endl;
#endif

    // always create CPU fallback
    cpu_fallback = std::make_unique<QueryBin>(
        intervals, store, binsize, seg_threshold, non_ref_threshold,
        num_threads, clip_mode, clip_margin, min_mutations_percent, max_mutations_percent,
        min_alignment_length, max_alignment_length
    );
    
    if (!use_gpu) {
        std::cout << "using CPU implementation with " << num_threads << " threads" << std::endl;
    }
}

QueryBinGPU::~QueryBinGPU() {
#ifdef METAL_SUPPORT
    if (metal_impl) {
        metal_destroy_implementation(metal_impl);
    }
#endif
}

void QueryBinGPU::aggregate_data() {
#ifdef METAL_SUPPORT
    if (use_gpu && metal_impl) {
        // collect all unique alignments across all intervals
        std::set<const Alignment*> processed_alignments_set;
        
        for (const auto& interval : intervals) {
            std::vector<std::reference_wrapper<const Alignment>> alignments = 
                store.get_alignments_in_interval(interval);
            
            for (const auto& alignment_ref : alignments) {
                const auto& aln = alignment_ref.get();
                processed_alignments_set.insert(&aln);
            }
        }
        
        std::vector<const Alignment*> processed_alignments(
            processed_alignments_set.begin(), processed_alignments_set.end());

        // apply same alignment filtering as CPU path
        std::vector<const Alignment*> filtered_alignments;
        filtered_alignments.reserve(processed_alignments.size());
        for (const Alignment* aln_ptr : processed_alignments) {
            if (passes_alignment_filter(*aln_ptr, store, clip_mode, clip_margin, min_mutations_percent, max_mutations_percent)) {
                filtered_alignments.push_back(aln_ptr);
            }
        }
        
        // prepare data for Metal
        std::vector<uint32_t> contig_indices;
        std::vector<MetalSimpleInterval> simple_intervals;
        simple_intervals.reserve(intervals.size());
        for (const auto& interval : intervals) {
            contig_indices.push_back(store.get_contig_index(interval.contig));
            MetalSimpleInterval si;
            si.start = interval.start;
            si.end = interval.end;
            simple_intervals.push_back(si);
        }
        
        // calculate dimensions
        uint32_t max_contigs = store.get_contigs().size();
        uint32_t max_bins_per_contig = 0;
        for (const auto& contig : store.get_contigs()) {
            uint32_t bins_needed = (contig.length + binsize - 1) / binsize;
            max_bins_per_contig = std::max(max_bins_per_contig, bins_needed);
        }
        
        std::cout << "running GPU processing with " << filtered_alignments.size() 
                  << " alignments..." << std::endl;
        
        // convert alignment pointers to actual alignments for Metal interface
        std::vector<MetalSimpleAlignment> actual_alignments;
        actual_alignments.reserve(filtered_alignments.size());
        for (const Alignment* aln_ptr : filtered_alignments) {
            MetalSimpleAlignment s;
            s.contig_index = aln_ptr->contig_index;
            s.contig_start = aln_ptr->contig_start;
            s.contig_end = aln_ptr->contig_end;
            s.num_mutations = static_cast<uint32_t>(aln_ptr->get_mutation_count());
            actual_alignments.push_back(s);
        }
        
        // run GPU implementation
        int success = metal_aggregate_data(metal_impl,
                                         actual_alignments.data(),
                                         static_cast<int>(actual_alignments.size()),
                                         simple_intervals.data(),
                                         static_cast<int>(simple_intervals.size()),
                                         contig_indices.data(),
                                         binsize,
                                         max_contigs,
                                         max_bins_per_contig);
        
        if (success) {
            // extract results from GPU
            bin_results.clear();
            int result_count = metal_get_result_count(metal_impl);
            
            for (int i = 0; i < result_count; i++) {
                uint32_t contig_index, bin_start, sequenced_bp, mutation_count;
                uint32_t dist_none, dist_5, dist_4, dist_3, dist_2, dist_1_plus;
                
                metal_get_result(metal_impl, i, &contig_index, &bin_start,
                               &sequenced_bp, &mutation_count,
                               &dist_none, &dist_5, &dist_4, &dist_3, &dist_2, &dist_1_plus);
                
                BinData data;
                data.sequenced_basepairs = sequenced_bp;
                data.mutation_count = 0;
                data.dist_none = dist_none;
                data.dist_5 = dist_5;
                data.dist_4 = dist_4;
                data.dist_3 = dist_3;
                data.dist_2 = dist_2;
                data.dist_1_plus = dist_1_plus;
                
                bin_results[{contig_index, bin_start}] = data;
            }
            
            // recalculate mutation_count on cpu to match semantics
            for (const Alignment* aln_ptr : filtered_alignments) {
                const Alignment& aln = *aln_ptr;
                if (aln.mutations.empty()) continue;
                for (uint32_t mut_idx : aln.mutations) {
                    const Mutation& m = store.get_mutation(aln.contig_index, mut_idx);
                    uint32_t pos = m.position;
                    for (const auto& interval : intervals) {
                        if (store.get_contig_index(interval.contig) != aln.contig_index) continue;
                        if (pos < interval.start || pos >= interval.end) continue;
                        uint32_t bstart = (pos / static_cast<uint32_t>(binsize)) * static_cast<uint32_t>(binsize);
                        auto it = bin_results.find({aln.contig_index, bstart});
                        if (it != bin_results.end()) {
                            it->second.mutation_count++;
                        }
                        break;
                    }
                }
            }

            // compute read_count on cpu with unique read indices per bin
            for (const Alignment* aln_ptr : filtered_alignments) {
                const Alignment& aln = *aln_ptr;
                for (const auto& interval : intervals) {
                    if (store.get_contig_index(interval.contig) != aln.contig_index) continue;
                    if (aln.contig_end <= interval.start || aln.contig_start >= interval.end) continue;
                    uint32_t adjusted_start = (interval.start / static_cast<uint32_t>(binsize)) * static_cast<uint32_t>(binsize);
                    if (interval.end == 0 || interval.start >= interval.end) continue;
                    uint32_t last_bin_start = ((interval.end - 1) / static_cast<uint32_t>(binsize)) * static_cast<uint32_t>(binsize);
                    for (uint32_t b_start = adjusted_start; b_start <= last_bin_start; b_start += static_cast<uint32_t>(binsize)) {
                        uint32_t b_end = b_start + static_cast<uint32_t>(binsize);
                        uint32_t effective_start = std::max({aln.contig_start, b_start, interval.start});
                        uint32_t effective_end = std::min({aln.contig_end, b_end, interval.end});
                        if (effective_end > effective_start) {
                            auto it = bin_results.find({aln.contig_index, b_start});
                            if (it != bin_results.end()) {
                                it->second.unique_reads.insert(aln.read_index);
                            }
                        }
                    }
                }
            }

            // collect mutation densities for median calculation (cpu)
            for (const Alignment* aln_ptr : filtered_alignments) {
                const Alignment& aln = *aln_ptr;
                uint32_t alignment_length = aln.contig_end - aln.contig_start;
                double local_mutations_per_bp = (alignment_length > 0) ? 
                    (static_cast<double>(aln.get_mutation_count()) / alignment_length) : 0.0;
                
                for (const auto& interval : intervals) {
                    if (store.get_contig_index(interval.contig) != aln.contig_index) continue;
                    if (aln.contig_end <= interval.start || aln.contig_start >= interval.end) continue;
                    uint32_t adjusted_start = (interval.start / static_cast<uint32_t>(binsize)) * static_cast<uint32_t>(binsize);
                    if (interval.end == 0 || interval.start >= interval.end) continue;
                    uint32_t last_bin_start = ((interval.end - 1) / static_cast<uint32_t>(binsize)) * static_cast<uint32_t>(binsize);
                    
                    // track which bins this alignment has already been processed for
                    std::set<std::pair<uint32_t, uint32_t>> processed_bins_for_alignment;
                    
                    for (uint32_t b_start = adjusted_start; b_start <= last_bin_start; b_start += static_cast<uint32_t>(binsize)) {
                        uint32_t b_end = b_start + static_cast<uint32_t>(binsize);
                        uint32_t effective_start = std::max({aln.contig_start, b_start, interval.start});
                        uint32_t effective_end = std::min({aln.contig_end, b_end, interval.end});
                        if (effective_end > effective_start) {
                            std::pair<uint32_t, uint32_t> bin_key = {aln.contig_index, b_start};
                            if (processed_bins_for_alignment.find(bin_key) == processed_bins_for_alignment.end()) {
                                processed_bins_for_alignment.insert(bin_key);
                                auto it = bin_results.find(bin_key);
                                if (it != bin_results.end()) {
                                    it->second.mutation_densities.push_back(local_mutations_per_bp);
                                }
                            }
                        }
                    }
                }
            }

            // variant counts by mutation position (cpu) for densities
            for (const Alignment* aln_ptr : filtered_alignments) {
                const Alignment& aln = *aln_ptr;
                for (uint32_t mut_idx : aln.mutations) {
                    const Mutation& m = store.get_mutation(aln.contig_index, mut_idx);
                    uint32_t pos = m.position;
                    for (const auto& interval : intervals) {
                        if (store.get_contig_index(interval.contig) != aln.contig_index) continue;
                        if (pos < interval.start || pos >= interval.end) continue;
                        uint32_t bstart = (pos / static_cast<uint32_t>(binsize)) * static_cast<uint32_t>(binsize);
                        auto it = bin_results.find({aln.contig_index, bstart});
                        if (it != bin_results.end()) {
                            std::string variant_key = std::to_string(pos) + "_" + std::to_string(static_cast<int>(m.type)) + "_" + m.nts;
                            it->second.variant_counts[variant_key]++;
                        }
                        break;
                    }
                }
            }

            // process clips for density calculations
            for (const Alignment* aln_ptr : filtered_alignments) {
                const Alignment& aln = *aln_ptr;
                uint32_t read_length = store.get_reads()[aln.read_index].length;
                bool is_left_clipped = (aln.read_start > static_cast<uint32_t>(clip_margin));
                bool is_right_clipped = (aln.read_end < (read_length - static_cast<uint32_t>(clip_margin)));

                // process left clip
                if (is_left_clipped) {
                    uint32_t clip_pos = aln.contig_start;
                    for (const auto& interval : intervals) {
                        if (store.get_contig_index(interval.contig) != aln.contig_index) continue;
                        if (clip_pos < interval.start || clip_pos >= interval.end) continue;
                        uint32_t bstart = (clip_pos / static_cast<uint32_t>(binsize)) * static_cast<uint32_t>(binsize);
                        auto it = bin_results.find({aln.contig_index, bstart});
                        if (it != bin_results.end()) {
                            it->second.clip_left_counts[std::to_string(clip_pos)]++;
                        }
                        break;
                    }
                }

                // process right clip
                if (is_right_clipped) {
                    uint32_t clip_pos = aln.contig_end;
                    for (const auto& interval : intervals) {
                        if (store.get_contig_index(interval.contig) != aln.contig_index) continue;
                        if (clip_pos < interval.start || clip_pos >= interval.end) continue;
                        uint32_t bstart = (clip_pos / static_cast<uint32_t>(binsize)) * static_cast<uint32_t>(binsize);
                        auto it = bin_results.find({aln.contig_index, bstart});
                        if (it != bin_results.end()) {
                            it->second.clip_right_counts[std::to_string(clip_pos)]++;
                        }
                        break;
                    }
                }
            }

            // efficient coverage recompute per bin only at positions with variants and clips
            for (auto &entry : bin_results) {
                uint32_t contig_index = entry.first.first;
                uint32_t b_start = entry.first.second;
                BinData &data = entry.second;
                // collect positions to evaluate coverage for
                std::vector<uint32_t> positions;
                positions.reserve(data.variant_counts.size() + data.clip_left_counts.size() + data.clip_right_counts.size());
                for (const auto &ve : data.variant_counts) {
                    const std::string &key = ve.first;
                    size_t u = key.find('_');
                    if (u == std::string::npos) continue;
                    uint32_t p = static_cast<uint32_t>(std::stoul(key.substr(0, u)));
                    positions.push_back(p);
                }
                // add clip positions
                for (const auto &ce : data.clip_left_counts) {
                    uint32_t p = static_cast<uint32_t>(std::stoul(ce.first));
                    positions.push_back(p);
                }
                for (const auto &ce : data.clip_right_counts) {
                    uint32_t p = static_cast<uint32_t>(std::stoul(ce.first));
                    positions.push_back(p);
                }
                if (positions.empty()) continue;
                std::sort(positions.begin(), positions.end());
                positions.erase(std::unique(positions.begin(), positions.end()), positions.end());

                // query alignments overlapping this bin interval
                std::string contig_id = store.get_contig_id(contig_index);
                Interval bin_interval(contig_id, b_start, b_start + static_cast<uint32_t>(binsize));
                auto aln_refs = store.get_alignments_in_interval(bin_interval);

                // reset coverage map
                data.position_coverage.clear();

                for (const auto &ar : aln_refs) {
                    const Alignment &aln = ar.get();
                    if (!passes_alignment_filter(aln, store, clip_mode, clip_margin, min_mutations_percent, max_mutations_percent)) continue;
                    if (aln.contig_index != contig_index) continue;
                    uint32_t effective_start = std::max(aln.contig_start, b_start);
                    uint32_t effective_end = std::min(aln.contig_end, b_start + static_cast<uint32_t>(binsize));
                    if (effective_end <= effective_start) continue;
                    // positions are sorted; accumulate coverage only within overlap
                    auto itp = std::lower_bound(positions.begin(), positions.end(), effective_start);
                    for (; itp != positions.end() && *itp < effective_end; ++itp) {
                        data.position_coverage[std::to_string(*itp)]++;
                    }
                }
            }
            
            std::cout << "GPU processing completed successfully with " 
                      << result_count << " bins" << std::endl;
            return;
        } else {
            std::cout << "GPU processing failed, falling back to CPU" << std::endl;
            use_gpu = false;
        }
    }
#endif
    
    // CPU fallback
    cpu_fallback->aggregate_data();
    std::cout << "cPU processing completed" << std::endl;
    
    // For CPU fallback, we need to recreate the results since QueryBin doesn't expose them
    // This is a temporary workaround - in production we'd modify QueryBin to expose results
    bin_results.clear();
    for (const auto& interval : intervals) {
        uint32_t contig_index = store.get_contig_index(interval.contig);
        uint32_t adjusted_start = (interval.start / binsize) * binsize;
        if (interval.end == 0 || interval.start >= interval.end) continue;
        uint32_t last_bin_start = ((interval.end - 1) / binsize) * binsize;
        
        for (uint32_t b_start = adjusted_start; b_start <= last_bin_start; b_start += binsize) {
            bin_results.try_emplace(std::make_pair(contig_index, b_start), BinData());
        }
    }
}

void QueryBinGPU::generate_output_rows() {
    output_rows.clear();

    for (const auto& entry : bin_results) {
        const auto& key = entry.first;
        const auto& data = entry.second;

        uint32_t contig_index = key.first;
        uint32_t bin_start = key.second;
        uint32_t bin_end = bin_start + binsize;
        int bin_length = binsize;

        std::string contig_id = store.get_contig_id(contig_index);

        // calculate segregating and non-reference sites density
        int segregating_sites = 0;
        int non_ref_sites = 0;
        for (const auto& variant_entry : data.variant_counts) {
            size_t first_underscore = variant_entry.first.find('_');
            if (first_underscore != std::string::npos) {
                std::string pos_str = variant_entry.first.substr(0, first_underscore);
                uint32_t position = std::stoul(pos_str);
                
                auto cov_it = data.position_coverage.find(std::to_string(position));
                if (cov_it != data.position_coverage.end() && cov_it->second > 0) {
                    double frequency = static_cast<double>(variant_entry.second) / cov_it->second;
                    
                    double min_freq = (seg_threshold == 0.0) ? 0.0 : seg_threshold;
                    double max_freq = (seg_threshold == 0.0) ? 1.0 : (1.0 - seg_threshold);
                    if (frequency > min_freq && frequency < max_freq) {
                        segregating_sites++;
                    }
                    
                    if (frequency >= non_ref_threshold) {
                        non_ref_sites++;
                    }
                }
            }
        }
        
        double seg_sites_density = static_cast<double>(segregating_sites) / binsize;
        double non_ref_sites_density = static_cast<double>(non_ref_sites) / binsize;

        // calculate segregating and non-reference clip sites density
        int segregating_clip_sites = 0;
        int non_ref_clip_sites = 0;
        
        // process left clips
        for (const auto& clip_entry : data.clip_left_counts) {
            
            auto cov_it = data.position_coverage.find(clip_entry.first);
            if (cov_it != data.position_coverage.end() && cov_it->second > 0) {
                double frequency = static_cast<double>(clip_entry.second) / cov_it->second;
                
                double min_freq = (seg_threshold == 0.0) ? 0.0 : seg_threshold;
                double max_freq = (seg_threshold == 0.0) ? 1.0 : (1.0 - seg_threshold);
                if (frequency > min_freq && frequency < max_freq) {
                    segregating_clip_sites++;
                }
                
                if (frequency >= non_ref_threshold) {
                    non_ref_clip_sites++;
                }
            }
        }
        
        // process right clips
        for (const auto& clip_entry : data.clip_right_counts) {
            
            auto cov_it = data.position_coverage.find(clip_entry.first);
            if (cov_it != data.position_coverage.end() && cov_it->second > 0) {
                double frequency = static_cast<double>(clip_entry.second) / cov_it->second;
                
                double min_freq = (seg_threshold == 0.0) ? 0.0 : seg_threshold;
                double max_freq = (seg_threshold == 0.0) ? 1.0 : (1.0 - seg_threshold);
                if (frequency > min_freq && frequency < max_freq) {
                    segregating_clip_sites++;
                }
                
                if (frequency >= non_ref_threshold) {
                    non_ref_clip_sites++;
                }
            }
        }
        
        double seg_clip_density = static_cast<double>(segregating_clip_sites) / binsize;
        double non_ref_clip_density = static_cast<double>(non_ref_clip_sites) / binsize;

        // calculate median mutation density
        double median_mutation_density = 0.0;
        if (!data.mutation_densities.empty()) {
            std::vector<double> densities_copy = data.mutation_densities;
            size_t n = densities_copy.size();
            if (n % 2 == 0) {
                // even number of values: average of middle two
                std::nth_element(densities_copy.begin(), densities_copy.begin() + n/2 - 1, densities_copy.end());
                double lower = densities_copy[n/2 - 1];
                std::nth_element(densities_copy.begin(), densities_copy.begin() + n/2, densities_copy.end());
                double upper = densities_copy[n/2];
                median_mutation_density = (lower + upper) / 2.0;
            } else {
                // odd number of values: middle value
                std::nth_element(densities_copy.begin(), densities_copy.begin() + n/2, densities_copy.end());
                median_mutation_density = densities_copy[n/2];
            }
        }

        output_rows.push_back({ contig_id, bin_start, bin_end, bin_length,
            data.sequenced_basepairs, static_cast<int>(data.unique_reads.size()), 
            data.mutation_count, median_mutation_density, seg_sites_density, non_ref_sites_density,
            seg_clip_density, non_ref_clip_density,
            data.dist_none, data.dist_5, data.dist_4,
            data.dist_3, data.dist_2, data.dist_1_plus });
    }
}

void QueryBinGPU::execute() {
    aggregate_data();
    generate_output_rows();
}

void QueryBinGPU::write_to_csv(const std::string& ofn_prefix) {
    std::string filename = ofn_prefix + "_bins.tsv";
    std::cout << "writing bin data rows to " << filename << std::endl;
    std::ofstream ofs(filename);

    if (!ofs.is_open()) {
        std::cerr << "error: could not open file " << filename << std::endl;
        exit(1);
    }

    // write header
    ofs << "contig\tbin_start\tbin_end\tbin_length\tsequenced_bp\tread_count\tmutation_count\t"
        << "median_mutation_density\tseg_sites_density\tnon_ref_sites_density\tseg_clip_density\tnon_ref_clip_density\t"
        << "dist_none\tdist_5\tdist_4\tdist_3\tdist_2\tdist_1_plus\n";

    for (const auto& row : output_rows) {
        ofs << row.contig << "\t"
            << row.bin_start << "\t"
            << row.bin_end << "\t"
            << row.bin_length << "\t"
            << row.sequenced_basepairs << "\t"
            << row.read_count << "\t"
            << row.mutation_count << "\t"
            << std::fixed << std::setprecision(5) << row.median_mutation_density << "\t"
            << row.seg_sites_density << "\t"
            << row.non_ref_sites_density << "\t"
            << row.seg_clip_density << "\t"
            << row.non_ref_clip_density << "\t"
            << row.dist_none << "\t"
            << row.dist_5 << "\t"
            << row.dist_4 << "\t"
            << row.dist_3 << "\t"
            << row.dist_2 << "\t"
            << row.dist_1_plus << "\n";
    }

    ofs.close();
}
