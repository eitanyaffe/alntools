#include "QueryBin.h"
#include <algorithm> // For std::min/max
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <thread>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

QueryBin::QueryBin(
    const std::vector<Interval>& intervals,
    const AlignmentStore& store,
    int binsize,
    double seg_threshold,
    double non_ref_threshold,
    int num_threads,
    ClipMode clip_mode,
    int clip_margin,
    double min_mutations_percent,
    double max_mutations_percent)
    : intervals(intervals)
    , store(store)
    , binsize(binsize)
    , seg_threshold(seg_threshold)
    , non_ref_threshold(non_ref_threshold)
    , num_threads(num_threads)
    , clip_mode(clip_mode)
    , clip_margin(clip_margin)
    , min_mutations_percent(min_mutations_percent)
    , max_mutations_percent(max_mutations_percent)
{
  if (binsize <= 0) {
    cerr << "error: binsize must be positive." << endl;
    exit(1);
  }
  // validate and correct seg_threshold
  if (seg_threshold < 0.0 || seg_threshold >= 0.5 || !std::isfinite(seg_threshold)) {
    cerr << "warning: invalid seg_threshold (" << seg_threshold << "), setting to default 0.2" << endl;
    this->seg_threshold = 0.2;
  }
  
  // validate and correct non_ref_threshold
  if (non_ref_threshold <= 0.5 || non_ref_threshold >= 1.0 || !std::isfinite(non_ref_threshold)) {
    cerr << "warning: invalid non_ref_threshold (" << non_ref_threshold << "), setting to default 0.9" << endl;
    this->non_ref_threshold = 0.9;
  }
  
  // set number of threads
  if (this->num_threads <= 0) {
    this->num_threads = std::thread::hardware_concurrency();
    if (this->num_threads <= 0) this->num_threads = 1;
  }
  
#ifdef _OPENMP
  omp_set_num_threads(this->num_threads);
#endif
}

void QueryBin::merge_bin_data(const std::map<std::pair<uint32_t, uint32_t>, BinData>& local_data)
{
  for (const auto& entry : local_data) {
    const auto& key = entry.first;
    const auto& local_bin = entry.second;
    
    auto it = bin_results.find(key);
    if (it != bin_results.end()) {
      // merge into existing bin
      BinData& global_bin = it->second;
      
      // sum simple counts
      global_bin.sequenced_basepairs += local_bin.sequenced_basepairs;
      global_bin.mutation_count += local_bin.mutation_count;
      global_bin.dist_none += local_bin.dist_none;
      global_bin.dist_5 += local_bin.dist_5;
      global_bin.dist_4 += local_bin.dist_4;
      global_bin.dist_3 += local_bin.dist_3;
      global_bin.dist_2 += local_bin.dist_2;
      global_bin.dist_1_plus += local_bin.dist_1_plus;
      
      // merge unique_reads set
      global_bin.unique_reads.insert(local_bin.unique_reads.begin(), local_bin.unique_reads.end());
      
      // merge variant_counts map
      for (const auto& variant_entry : local_bin.variant_counts) {
        global_bin.variant_counts[variant_entry.first] += variant_entry.second;
      }
      
      // merge position_coverage map
      for (const auto& pos_entry : local_bin.position_coverage) {
        global_bin.position_coverage[pos_entry.first] += pos_entry.second;
      }
      
      // merge clip_left_counts map
      for (const auto& clip_entry : local_bin.clip_left_counts) {
        global_bin.clip_left_counts[clip_entry.first] += clip_entry.second;
      }
      
      // merge clip_right_counts map
      for (const auto& clip_entry : local_bin.clip_right_counts) {
        global_bin.clip_right_counts[clip_entry.first] += clip_entry.second;
      }
    } else {
      // create new bin entry
      bin_results[key] = local_bin;
    }
  }
}

void QueryBin::process_single_alignment(const Alignment& aln, std::map<std::pair<uint32_t, uint32_t>, BinData>& target_bin_results)
{
  // apply alignment filtering
  if (!passes_alignment_filter(aln, store, clip_mode, clip_margin, min_mutations_percent, max_mutations_percent)) {
    return;
  }
  
  // calculate alignment length for mutation distance calculation
  uint32_t alignment_length = aln.contig_end - aln.contig_start;
  double local_mutations_per_bp = (alignment_length > 0) ? 
    (static_cast<double>(aln.mutations.size()) / alignment_length) : 0.0;

  // process sequenced basepairs and mutation rate categories for this alignment across all relevant intervals
  for (const auto& interval : intervals) {
    if (aln.contig_index != store.get_contig_index(interval.contig)) {
      continue;
    }

    // skip if alignment doesn't overlap this interval
    if (aln.contig_end <= interval.start || aln.contig_start >= interval.end) {
      continue;
    }

    uint32_t adjusted_start = (interval.start / binsize) * binsize;
    uint32_t last_bin_start = ((interval.end - 1) / binsize) * binsize;

    // check overlap with each bin in this interval
    for (uint32_t b_start = adjusted_start; b_start <= last_bin_start; b_start += binsize) {
      uint32_t b_end = b_start + binsize;

      // calculate overlap considering alignment, bin, AND original interval boundaries
      uint32_t effective_start = std::max({ aln.contig_start, b_start, interval.start });
      uint32_t effective_end = std::min({ aln.contig_end, b_end, interval.end });

      int overlap_length = (effective_end > effective_start) ? (effective_end - effective_start) : 0;

      if (overlap_length > 0) {
        auto it = target_bin_results.find({ aln.contig_index, b_start });
        if (it != target_bin_results.end()) {
          it->second.sequenced_basepairs += overlap_length;
          it->second.unique_reads.insert(aln.read_index); // track unique reads
          
          // track coverage for segregating sites calculation
          for (uint32_t pos = effective_start; pos < effective_end; pos++) {
            it->second.position_coverage[std::to_string(pos)]++;
          }
        }
      }
    }
  }

  // process mutations for this alignment only once
  for (uint32_t mutation_index : aln.mutations) {
    // fetch the mutation object
    const Mutation& mutation = store.get_mutation(aln.contig_index, mutation_index);
    uint32_t mutation_contig_pos = mutation.position;
    uint32_t mutation_bin_start = (mutation_contig_pos / binsize) * binsize;

    // check if this mutation falls within any of our intervals
    for (const auto& interval : intervals) {
      if (aln.contig_index != store.get_contig_index(interval.contig)) {
        continue;
      }

      // check if mutation is within this interval
      if (mutation_contig_pos >= interval.start && mutation_contig_pos < interval.end) {
        auto it = target_bin_results.find({ aln.contig_index, mutation_bin_start });
        if (it != target_bin_results.end()) {
          it->second.mutation_count++;
          
          // create variant key for segregating sites
          std::string variant_key = std::to_string(mutation_contig_pos) + "_" + 
            std::to_string(static_cast<int>(mutation.type)) + "_" + mutation.nts;
          it->second.variant_counts[variant_key]++;
          
          break; // only count the mutation once even if it's in multiple overlapping intervals
        }
      }
    }
  }

  // process clips for this alignment
  uint32_t read_length = store.get_reads()[aln.read_index].length;
  bool is_left_clipped = (aln.read_start > static_cast<uint32_t>(clip_margin));
  bool is_right_clipped = (aln.read_end < (read_length - static_cast<uint32_t>(clip_margin)));

  // process left clip
  if (is_left_clipped) {
    uint32_t clip_contig_pos = aln.contig_start;
    uint32_t clip_bin_start = (clip_contig_pos / binsize) * binsize;
    
    for (const auto& interval : intervals) {
      if (aln.contig_index != store.get_contig_index(interval.contig)) {
        continue;
      }
      
      // check if clip position is within this interval
      if (clip_contig_pos >= interval.start && clip_contig_pos < interval.end) {
        auto it = target_bin_results.find({ aln.contig_index, clip_bin_start });
        if (it != target_bin_results.end()) {
          it->second.clip_left_counts[std::to_string(clip_contig_pos)]++;
          break;
        }
      }
    }
  }

  // process right clip
  if (is_right_clipped) {
    uint32_t clip_contig_pos = aln.contig_end;
    uint32_t clip_bin_start = (clip_contig_pos / binsize) * binsize;
    
    for (const auto& interval : intervals) {
      if (aln.contig_index != store.get_contig_index(interval.contig)) {
        continue;
      }
      
      // check if clip position is within this interval
      if (clip_contig_pos >= interval.start && clip_contig_pos < interval.end) {
        auto it = target_bin_results.find({ aln.contig_index, clip_bin_start });
        if (it != target_bin_results.end()) {
          it->second.clip_right_counts[std::to_string(clip_contig_pos)]++;
          break;
        }
      }
    }
  }

  // categorize this alignment by mutation rate for all bins it overlaps
  std::set<std::pair<uint32_t, uint32_t>> processed_bins_for_alignment;
  for (const auto& interval : intervals) {
    if (aln.contig_index != store.get_contig_index(interval.contig)) {
      continue;
    }

    // skip if alignment doesn't overlap this interval
    if (aln.contig_end <= interval.start || aln.contig_start >= interval.end) {
      continue;
    }

    uint32_t adjusted_start = (interval.start / binsize) * binsize;
    uint32_t last_bin_start = ((interval.end - 1) / binsize) * binsize;

    // check overlap with each bin in this interval
    for (uint32_t b_start = adjusted_start; b_start <= last_bin_start; b_start += binsize) {
      uint32_t b_end = b_start + binsize;

      // calculate overlap considering alignment, bin, AND original interval boundaries
      uint32_t effective_start = std::max({ aln.contig_start, b_start, interval.start });
      uint32_t effective_end = std::min({ aln.contig_end, b_end, interval.end });

      int overlap_length = (effective_end > effective_start) ? (effective_end - effective_start) : 0;

      if (overlap_length > 0) {
        // only process each bin once per alignment
        std::pair<uint32_t, uint32_t> bin_key = { aln.contig_index, b_start };
        if (processed_bins_for_alignment.find(bin_key) == processed_bins_for_alignment.end()) {
          processed_bins_for_alignment.insert(bin_key);
          
          auto it = target_bin_results.find(bin_key);
          if (it != target_bin_results.end()) {
            // categorize alignment by mutation distance (per bp)
            if (aln.mutations.size() == 0) {
              it->second.dist_none++;
            } else if (local_mutations_per_bp < 1e-4) {
              it->second.dist_5++;  // 1e-5 to 1e-4 per bp
            } else if (local_mutations_per_bp < 1e-3) {
              it->second.dist_4++;  // 1e-4 to 1e-3 per bp
            } else if (local_mutations_per_bp < 1e-2) {
              it->second.dist_3++;  // 1e-3 to 1e-2 per bp
            } else if (local_mutations_per_bp < 1e-1) {
              it->second.dist_2++;  // 1e-2 to 1e-1 per bp
            } else {
              it->second.dist_1_plus++;  // above 1e-1 per bp
            }
          }
        }
      }
    }
  }
}

void QueryBin::aggregate_data()
{
  bin_results.clear();

  // First pass: initialize all relevant bins across all intervals
  for (const auto& interval : intervals) {
    uint32_t contig_index = store.get_contig_index(interval.contig);

    // Calculate relevant bin range based on original interval
    uint32_t adjusted_start = (interval.start / binsize) * binsize;
    // Handle edge case where end is 0 or interval is empty
    if (interval.end == 0 || interval.start >= interval.end)
      continue;
    uint32_t last_bin_start = ((interval.end - 1) / binsize) * binsize;

    // Initialize relevant bins in the map
    for (uint32_t b_start = adjusted_start; b_start <= last_bin_start; b_start += binsize) {
      bin_results.try_emplace({ contig_index, b_start }, BinData());
    }
  }

  // Collect all unique alignments across all intervals
  std::set<const Alignment*> processed_alignments_set;
  
  for (const auto& interval : intervals) {
    // Get alignments overlapping this interval
    std::vector<std::reference_wrapper<const Alignment>> alignments = store.get_alignments_in_interval(interval);
    
    for (const auto& alignment_ref : alignments) {
      const auto& aln = alignment_ref.get();
      processed_alignments_set.insert(&aln);
    }
  }

  // Convert to vector for indexed access in OpenMP
  std::vector<const Alignment*> processed_alignments(processed_alignments_set.begin(), processed_alignments_set.end());

  // Start timing
  auto start_time = std::chrono::high_resolution_clock::now();

  // Second pass: process each alignment only once (parallelized)
#ifdef _OPENMP
  #pragma omp parallel
  {
    #pragma omp single
    {
      cout << "using " << omp_get_num_threads() << " threads for " << processed_alignments.size() << " alignments" << endl;
    }
    // thread-local bin results
    std::map<std::pair<uint32_t, uint32_t>, BinData> local_bin_results;
    
    // initialize local bins with empty data
    for (const auto& entry : bin_results) {
      local_bin_results[entry.first] = BinData();
    }
    
    #pragma omp for
    for (size_t i = 0; i < processed_alignments.size(); ++i) {
      const Alignment* aln_ptr = processed_alignments[i];
      process_single_alignment(*aln_ptr, local_bin_results);
    }
    
    // merge local results into global bin_results
    #pragma omp critical
    {
      merge_bin_data(local_bin_results);
    }
  }
#else
  cout << "processing " << processed_alignments.size() << " alignments sequentially" << endl;
  for (const Alignment* aln_ptr : processed_alignments) {
    process_single_alignment(*aln_ptr, bin_results);
  }
#endif

  // End timing and report
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
  cout << "binning completed in " << duration.count() << " ms" << endl;
}

void QueryBin::generate_output_rows()
{
  output_rows.clear();

  for (const auto& entry : bin_results) {
    const auto& key = entry.first;
    const auto& data = entry.second;

    uint32_t contig_index = key.first;
    uint32_t bin_start = key.second;
    uint32_t bin_end = bin_start + binsize; // Bin end is standard
    int bin_length = binsize; // Bin length is standard

    string contig_id = store.get_contig_id(contig_index);

    // calculate segregating and non-reference sites density
    int segregating_sites = 0;
    int non_ref_sites = 0;
    for (const auto& variant_entry : data.variant_counts) {
      // extract position from variant key
      size_t first_underscore = variant_entry.first.find('_');
      if (first_underscore != string::npos) {
        string pos_str = variant_entry.first.substr(0, first_underscore);
        uint32_t position = std::stoul(pos_str);
        
        // get coverage at this position
        auto cov_it = data.position_coverage.find(std::to_string(position));
        if (cov_it != data.position_coverage.end() && cov_it->second > 0) {
          double frequency = static_cast<double>(variant_entry.second) / cov_it->second;
          
          // check if segregating
          double min_freq = (seg_threshold == 0.0) ? 0.0 : seg_threshold;
          double max_freq = (seg_threshold == 0.0) ? 1.0 : (1.0 - seg_threshold);
          if (frequency > min_freq && frequency < max_freq) {
            segregating_sites++;
          }
          
          // check if non-reference (high frequency)
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
      
      // get coverage at this position
      auto cov_it = data.position_coverage.find(clip_entry.first);
      if (cov_it != data.position_coverage.end() && cov_it->second > 0) {
        double frequency = static_cast<double>(clip_entry.second) / cov_it->second;
        
        // check if segregating
        double min_freq = (seg_threshold == 0.0) ? 0.0 : seg_threshold;
        double max_freq = (seg_threshold == 0.0) ? 1.0 : (1.0 - seg_threshold);
        if (frequency > min_freq && frequency < max_freq) {
          segregating_clip_sites++;
        }
        
        // check if non-reference (high frequency)
        if (frequency >= non_ref_threshold) {
          non_ref_clip_sites++;
        }
      }
    }
    
    // process right clips
    for (const auto& clip_entry : data.clip_right_counts) {
      
      // get coverage at this position
      auto cov_it = data.position_coverage.find(clip_entry.first);
      if (cov_it != data.position_coverage.end() && cov_it->second > 0) {
        double frequency = static_cast<double>(clip_entry.second) / cov_it->second;
        
        // check if segregating
        double min_freq = (seg_threshold == 0.0) ? 0.0 : seg_threshold;
        double max_freq = (seg_threshold == 0.0) ? 1.0 : (1.0 - seg_threshold);
        if (frequency > min_freq && frequency < max_freq) {
          segregating_clip_sites++;
        }
        
        // check if non-reference (high frequency)
        if (frequency >= non_ref_threshold) {
          non_ref_clip_sites++;
        }
      }
    }
    
    double seg_clip_density = static_cast<double>(segregating_clip_sites) / binsize;
    double non_ref_clip_density = static_cast<double>(non_ref_clip_sites) / binsize;

    output_rows.push_back({ contig_id, bin_start, bin_end, bin_length,
        data.sequenced_basepairs, static_cast<int>(data.unique_reads.size()), data.mutation_count, seg_sites_density, non_ref_sites_density,
        seg_clip_density, non_ref_clip_density,
        data.dist_none, data.dist_5, data.dist_4,
        data.dist_3, data.dist_2, data.dist_1_plus });
  }
}

// function to write the generated rows to a file
void QueryBin::write_to_csv(const std::string& ofn_prefix)
{
  string filename = ofn_prefix + "_bins.tsv";
  cout << "writing bin data rows to " << filename << endl;
  ofstream ofs(filename);

  if (!ofs.is_open()) {
    // Consider throwing an exception instead of exiting
    cerr << "error: could not open file " << filename << endl;
    exit(1);
  }

  // Write header with columns
  ofs << "contig\tbin_start\tbin_end\tbin_length\tsequenced_bp\tread_count\tmutation_count\t"
      << "seg_sites_density\tnon_ref_sites_density\tseg_clip_density\tnon_ref_clip_density\t"
      << "dist_none\tdist_5\tdist_4\tdist_3\tdist_2\tdist_1_plus\n";

  for (const auto& row : output_rows) {
    ofs << row.contig << "\t"
        << row.bin_start << "\t"
        << row.bin_end << "\t"
        << row.bin_length << "\t"
        << row.sequenced_basepairs << "\t"
        << row.read_count << "\t"
        << row.mutation_count << "\t"
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

void QueryBin::execute()
{
  aggregate_data();
  generate_output_rows();
}
