#include "SegmentationManager.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>

using namespace std;

SegmentationManager::SegmentationManager(const map<string, AlignmentStore>& stores,
                                       int max_margin,
                                       int min_anchor_length,
                                       int min_dangle_length,
                                       double max_anchor_mutations_percent,
                                       int min_alignment_distance,
                                       int min_breakpoint_read_support,
                                       double min_breakpoint_frequency,
                                       int min_segment_length)
    : stores(stores)
    , max_margin(max_margin)
    , min_anchor_length(min_anchor_length)
    , min_dangle_length(min_dangle_length)
    , max_anchor_mutations_percent(max_anchor_mutations_percent)
    , min_alignment_distance(min_alignment_distance)
    , min_breakpoint_read_support(min_breakpoint_read_support)
    , min_breakpoint_frequency(min_breakpoint_frequency)
    , min_segment_length(min_segment_length)
{
    for (const auto& entry : stores) {
        library_ids.push_back(entry.first);
    }
}

void SegmentationManager::execute()
{
    cout << "executing segmentation analysis across " << library_ids.size() << " libraries..." << endl;
    
    detect_breakpoints_per_library();
    aggregate_breakpoints_step();
    calculate_coverage();
    filter_breakpoints();
    generate_segments();
    
    cout << "segmentation analysis completed" << endl;
}

void SegmentationManager::detect_breakpoints_per_library()
{
    cout << "detecting breakpoints per library..." << endl;
    
    for (const string& lib_id : library_ids) {
        cout << "--------------------------------" << endl;
        cout << "processing library: " << lib_id << endl;
        
        Segmentation detector(stores.at(lib_id), max_margin, min_anchor_length,
                            min_dangle_length, max_anchor_mutations_percent,
                            min_alignment_distance);
        
        vector<ReadBreakpoint>& breakpoints = lib_breakpoints[lib_id];
        detector.detect_breakpoints(lib_id, breakpoints);
        
        cout << "library " << lib_id << ": found " << breakpoints.size() << " breakpoints" << endl;
    }
}

void SegmentationManager::aggregate_breakpoints_step()
{
    cout << "aggregating breakpoints across libraries..." << endl;
    
    // collect all breakpoints from all libraries
    vector<ReadBreakpoint> all_breakpoints;
    for (const auto& entry : lib_breakpoints) {
        const vector<ReadBreakpoint>& lib_bps = entry.second;
        
        for (const ReadBreakpoint& bp : lib_bps) {
            all_breakpoints.push_back(bp);
        }
    }
    
    if (all_breakpoints.empty()) {
        cout << "no breakpoints to aggregate" << endl;
        return;
    }
    
    cout << "total breakpoints before aggregation: " << all_breakpoints.size() << endl;
    
    // sort by contig, type, coord
    sort(all_breakpoints.begin(), all_breakpoints.end(),
         [](const ReadBreakpoint& a, const ReadBreakpoint& b) {
             if (a.contig_id != b.contig_id) return a.contig_id < b.contig_id;
             if (a.type != b.type) return a.type < b.type;
             return a.coord < b.coord;
         });
    
    // cluster breakpoints within margin
    aggregate_breakpoints.clear();
    int breakpoint_counter = 1;
    
    size_t i = 0;
    while (i < all_breakpoints.size()) {
        const ReadBreakpoint& first_bp = all_breakpoints[i];
        
        // collect all breakpoints within margin
        vector<ReadBreakpoint> cluster;
        cluster.push_back(first_bp);
        
        size_t j = i + 1;
        while (j < all_breakpoints.size()) {
            const ReadBreakpoint& next_bp = all_breakpoints[j];
            
            if (next_bp.contig_id != first_bp.contig_id || 
                next_bp.type != first_bp.type) {
                break;
            }
            
            int32_t coord_diff = abs(static_cast<int32_t>(next_bp.coord) - 
                                    static_cast<int32_t>(first_bp.coord));
            if (coord_diff > max_margin) {
                break;
            }
            
            cluster.push_back(next_bp);
            j++;
        }
        
        // find coord with most support (mode of coords)
        map<uint32_t, int> coord_counts;
        for (const ReadBreakpoint& bp : cluster) {
            coord_counts[bp.coord]++;
        }
        
        uint32_t best_coord = first_bp.coord;
        int max_count = 0;
        for (const auto& entry : coord_counts) {
            if (entry.second > max_count) {
                max_count = entry.second;
                best_coord = entry.first;
            }
        }
        
        // create aggregate breakpoint
        string bp_id = "b" + to_string(breakpoint_counter++);
        AggregateBreakpoint agg_bp(bp_id, first_bp.contig_id, best_coord, first_bp.type);
        agg_bp.read_support = cluster.size();
        
        // count support per library
        for (const ReadBreakpoint& bp : cluster) {
            // find which library this breakpoint came from
            for (const auto& lib_entry : lib_breakpoints) {
                const string& lib_id = lib_entry.first;
                const vector<ReadBreakpoint>& lib_bps = lib_entry.second;
                
                for (const ReadBreakpoint& lib_bp : lib_bps) {
                    if (lib_bp.read_id == bp.read_id && 
                        lib_bp.contig_id == bp.contig_id &&
                        lib_bp.coord == bp.coord &&
                        lib_bp.type == bp.type) {
                        agg_bp.lib_support[lib_id]++;
                        break;
                    }
                }
            }
        }
        
        aggregate_breakpoints.push_back(agg_bp);
        i = j;
    }
    
    cout << "created " << aggregate_breakpoints.size() << " aggregated breakpoints" << endl;
}

void SegmentationManager::calculate_coverage()
{
    cout << "calculating coverage for breakpoints..." << endl;
    
    for (AggregateBreakpoint& bp : aggregate_breakpoints) {
        int total_coverage = 0;
        
        for (const string& lib_id : library_ids) {
            int lib_cov = calculate_breakpoint_coverage(bp, lib_id);
            bp.lib_coverage[lib_id] = lib_cov;
            total_coverage += lib_cov;
        }
        
        // calculate frequency
        if (total_coverage > 0) {
            bp.frequency = static_cast<double>(bp.read_support) / static_cast<double>(total_coverage);
        } else {
            bp.frequency = 0.0;
        }
    }
    
    cout << "coverage calculation completed" << endl;
}

int SegmentationManager::calculate_breakpoint_coverage(const AggregateBreakpoint& bp, 
                                                      const string& lib_id) const
{
    const AlignmentStore& store = stores.at(lib_id);
    
    if (!store.has_contig_id(bp.contig_id)) {
        return 0;
    }
    
    // create interval around breakpoint coordinate
    uint32_t coord_0based = bp.coord > 0 ? bp.coord - 1 : 0;
    Interval interval(bp.contig_id, coord_0based, coord_0based + 1);
    
    auto alignments = store.get_alignments_intersecting_interval(interval);
    
    // collect unique read IDs that pass mutation filter
    set<string> reads;
    
    for (const auto& alignment_ref : alignments) {
        const Alignment& aln = alignment_ref.get();
        
        uint32_t alignment_length = aln.read_end - aln.read_start;
        double mutation_percent = (static_cast<double>(aln.get_mutation_count()) / 
            static_cast<double>(alignment_length)) * 100.0;
        
        if (mutation_percent <= max_anchor_mutations_percent) {
            string read_id = store.get_read_id(aln.read_index);
            reads.insert(read_id);
        }
    }
    
    return static_cast<int>(reads.size());
}

void SegmentationManager::filter_breakpoints()
{
    cout << "filtering breakpoints..." << endl;
    
    int selected_count = 0;
    int selected_global = 0;
    int selected_per_lib = 0;
    
    for (AggregateBreakpoint& bp : aggregate_breakpoints) {
        bool passes_global = (bp.read_support >= min_breakpoint_read_support && 
                             bp.frequency >= min_breakpoint_frequency);
        
        // check per-library thresholds
        bool passes_per_lib = false;
        for (const string& lib_id : library_ids) {
            auto support_it = bp.lib_support.find(lib_id);
            auto coverage_it = bp.lib_coverage.find(lib_id);
            
            if (support_it == bp.lib_support.end() || coverage_it == bp.lib_coverage.end()) {
                continue;
            }
            
            int lib_support = support_it->second;
            int lib_coverage = coverage_it->second;
            
            if (lib_coverage == 0) {
                continue;
            }
            
            double lib_frequency = static_cast<double>(lib_support) / static_cast<double>(lib_coverage);
            
            if (lib_support >= min_breakpoint_read_support && 
                lib_frequency >= min_breakpoint_frequency) {
                passes_per_lib = true;
                break;
            }
        }
        
        if (passes_global || passes_per_lib) {
            bp.selected = true;
            selected_count++;
            if (passes_global) selected_global++;
            if (passes_per_lib) selected_per_lib++;
        } else {
            bp.selected = false;
        }
    }
    
    cout << "selected " << selected_count << " / " << aggregate_breakpoints.size() 
         << " breakpoints" << endl;
    cout << "  (global threshold: " << selected_global << ", per-library threshold: " 
         << selected_per_lib << ")" << endl;
}

void SegmentationManager::generate_segments()
{
    cout << "generating segments..." << endl;
    
    segments.clear();
    
    // get all selected breakpoints
    vector<AggregateBreakpoint*> selected_bps;
    for (AggregateBreakpoint& bp : aggregate_breakpoints) {
        if (bp.selected) {
            selected_bps.push_back(&bp);
        }
    }
    
    if (selected_bps.empty()) {
        cout << "no selected breakpoints - no segments to generate" << endl;
        return;
    }
    
    // group by contig
    map<string, vector<AggregateBreakpoint*>> contig_breakpoints;
    for (AggregateBreakpoint* bp : selected_bps) {
        contig_breakpoints[bp->contig_id].push_back(bp);
    }
    
    int segment_counter = 1;
    
    for (auto& entry : contig_breakpoints) {
        const string& contig_id = entry.first;
        vector<AggregateBreakpoint*>& bps = entry.second;
        
        // sort by coordinate
        sort(bps.begin(), bps.end(),
             [](const AggregateBreakpoint* a, const AggregateBreakpoint* b) {
                 return a->coord < b->coord;
             });
        
        // get contig length from first store that has it
        uint32_t contig_length = 0;
        for (const auto& store_entry : stores) {
            const AlignmentStore& store = store_entry.second;
            if (store.has_contig_id(contig_id)) {
                size_t contig_idx = store.get_contig_index(contig_id);
                contig_length = store.get_contigs()[contig_idx].length;
                break;
            }
        }
        
        if (contig_length == 0) {
            continue;
        }
        
        // filter breakpoints based on min_segment_length
        vector<AggregateBreakpoint*> filtered_bps;
        uint32_t last_accepted_coord = 1;
        
        for (AggregateBreakpoint* bp : bps) {
            uint32_t segment_length = bp->coord - last_accepted_coord;
            
            if (segment_length >= static_cast<uint32_t>(min_segment_length)) {
                // check that this breakpoint also leaves enough space to contig end
                uint32_t remaining_length = contig_length - bp->coord;
                
                if (remaining_length >= static_cast<uint32_t>(min_segment_length)) {
                    filtered_bps.push_back(bp);
                    bp->is_segment_break = true;
                    last_accepted_coord = bp->coord;
                }
            }
        }
        
        // create segment from start to first breakpoint (if any)
        if (!filtered_bps.empty() && filtered_bps[0]->coord > 1) {
            string seg_id = "s" + to_string(segment_counter++);
            Segment seg(seg_id, contig_id, 1, filtered_bps[0]->coord);
            segments.push_back(seg);
        } else if (filtered_bps.empty()) {
            // no breakpoints passed filter, create single segment for entire contig
            string seg_id = "s" + to_string(segment_counter++);
            Segment seg(seg_id, contig_id, 1, contig_length);
            segments.push_back(seg);
        }
        
        // create segments between filtered breakpoints
        for (size_t i = 0; i + 1 < filtered_bps.size(); ++i) {
            if (filtered_bps[i]->coord == filtered_bps[i + 1]->coord) {
                continue;
            }
            string seg_id = "s" + to_string(segment_counter++);
            uint32_t seg_start = filtered_bps[i]->coord + 1;
            uint32_t seg_end = filtered_bps[i + 1]->coord;
            Segment seg(seg_id, contig_id, seg_start, seg_end);
            segments.push_back(seg);
        }
        
        // create segment from last breakpoint to end
        if (!filtered_bps.empty() && filtered_bps.back()->coord < contig_length) {
            string seg_id = "s" + to_string(segment_counter++);
            uint32_t seg_start = filtered_bps.back()->coord + 1;
            Segment seg(seg_id, contig_id, seg_start, contig_length);
            segments.push_back(seg);
        }
    }
    
    // report segments shorter than threshold
    int short_segments = 0;
    uint32_t min_length = UINT32_MAX;
    
    for (const Segment& seg : segments) {
        if (seg.length < static_cast<uint32_t>(min_segment_length)) {
            short_segments++;
            if (seg.length < min_length) {
                min_length = seg.length;
            }
        }
    }
    
    cout << "generated " << segments.size() << " segments" << endl;
    
    if (short_segments > 0) {
        cout << "warning: " << short_segments << " segments are shorter than min_segment_length (" 
             << min_segment_length << " bp), minimum length: " << min_length << " bp" << endl;
    }
}

void SegmentationManager::write_to_csv(const string& odir)
{
    cout << "writing output files to directory: " << odir << endl;
    
    write_read_breakpoints_file(odir);
    write_aggregate_breakpoints_file(odir);
    write_support_matrix_file(odir);
    write_coverage_matrix_file(odir);
    write_segments_file(odir);
    
    cout << "output files written" << endl;
}

void SegmentationManager::write_read_breakpoints_file(const string& odir)
{
    string filename = odir + "/read_breakpoints.txt";
    ofstream file;
    safe_open_file_for_writing(filename, file);
    
    // header
    file << "lib_id\tread_id\tcontig\tcoord\ttype\tanchor_length\tanchor_mutations\tdangle_length" << endl;
    
    // write all read breakpoints
    for (const auto& entry : lib_breakpoints) {
        const string& lib_id = entry.first;
        const vector<ReadBreakpoint>& breakpoints = entry.second;
        
        for (const ReadBreakpoint& bp : breakpoints) {
            file << lib_id << "\t"
                 << bp.read_id << "\t"
                 << bp.contig_id << "\t"
                 << bp.coord << "\t"
                 << bp.type << "\t"
                 << bp.anchor_length << "\t"
                 << bp.anchor_mutations << "\t"
                 << bp.dangle_length << endl;
        }
    }
    
    file.close();
    cout << "wrote " << filename << endl;
}

void SegmentationManager::write_aggregate_breakpoints_file(const string& odir)
{
    string filename = odir + "/breakpoints.txt";
    ofstream file;
    safe_open_file_for_writing(filename, file);
    
    // header
    file << "breakpoint_id\tcontig\tcoord\ttype\tread_support\tfrequency\tselected\tis_segment_break" << endl;
    
    for (const AggregateBreakpoint& bp : aggregate_breakpoints) {
        file << bp.breakpoint_id << "\t"
             << bp.contig_id << "\t"
             << bp.coord << "\t"
             << bp.type << "\t"
             << bp.read_support << "\t"
             << fixed << setprecision(4) << bp.frequency << "\t"
             << (bp.selected ? "T" : "F") << "\t"
             << (bp.is_segment_break ? "T" : "F") << endl;
    }
    
    file.close();
    cout << "wrote " << filename << endl;
}

void SegmentationManager::write_support_matrix_file(const string& odir)
{
    string filename = odir + "/breakpoints_support.txt";
    ofstream file;
    safe_open_file_for_writing(filename, file);
    
    // header
    file << "breakpoint_id";
    for (const string& lib_id : library_ids) {
        file << "\t" << lib_id;
    }
    file << endl;
    
    // data rows
    for (const AggregateBreakpoint& bp : aggregate_breakpoints) {
        file << bp.breakpoint_id;
        for (const string& lib_id : library_ids) {
            auto it = bp.lib_support.find(lib_id);
            int count = (it != bp.lib_support.end()) ? it->second : 0;
            file << "\t" << count;
        }
        file << endl;
    }
    
    file.close();
    cout << "wrote " << filename << endl;
}

void SegmentationManager::write_coverage_matrix_file(const string& odir)
{
    string filename = odir + "/breakpoints_coverage.txt";
    ofstream file;
    safe_open_file_for_writing(filename, file);
    
    // header
    file << "breakpoint_id";
    for (const string& lib_id : library_ids) {
        file << "\t" << lib_id;
    }
    file << endl;
    
    // data rows
    for (const AggregateBreakpoint& bp : aggregate_breakpoints) {
        file << bp.breakpoint_id;
        for (const string& lib_id : library_ids) {
            auto it = bp.lib_coverage.find(lib_id);
            int count = (it != bp.lib_coverage.end()) ? it->second : 0;
            file << "\t" << count;
        }
        file << endl;
    }
    
    file.close();
    cout << "wrote " << filename << endl;
}

void SegmentationManager::write_segments_file(const string& odir)
{
    string filename = odir + "/segments.txt";
    ofstream file;
    safe_open_file_for_writing(filename, file);
    
    // header
    file << "segment\tcontig\tstart\tend\tlength" << endl;
    
    for (const Segment& seg : segments) {
        file << seg.segment_id << "\t"
             << seg.contig_id << "\t"
             << seg.start << "\t"
             << seg.end << "\t"
             << seg.length << endl;
    }
    
    file.close();
    cout << "wrote " << filename << endl;
}

