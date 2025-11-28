#include "SegmentationManager.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>

using namespace std;

SegmentationManager::SegmentationManager(const map<string, AlignmentStore>& stores,
                                       const string& ifn_contig_table,
                                       int max_margin,
                                       int min_anchor_length,
                                       int min_dangle_length,
                                       double max_anchor_mutations_percent,
                                       int min_alignment_distance,
                                       int min_breakpoint_read_support,
                                       double min_breakpoint_frequency,
                                       int min_segment_length,
                                       int max_segment_length)
    : stores(stores)
    , ifn_contig_table(ifn_contig_table)
    , max_margin(max_margin)
    , min_anchor_length(min_anchor_length)
    , min_dangle_length(min_dangle_length)
    , max_anchor_mutations_percent(max_anchor_mutations_percent)
    , min_alignment_distance(min_alignment_distance)
    , min_breakpoint_read_support(min_breakpoint_read_support)
    , min_breakpoint_frequency(min_breakpoint_frequency)
    , min_segment_length(min_segment_length)
    , max_segment_length(max_segment_length)
{
    for (const auto& entry : stores) {
        library_ids.push_back(entry.first);
    }
}

void SegmentationManager::execute()
{
    cout << "executing segmentation analysis across " << library_ids.size() << " libraries..." << endl;
    
    load_contig_table();
    detect_breakpoints_per_library();
    aggregate_breakpoints_step();
    calculate_coverage();
    filter_breakpoints();
    generate_segments();
    
    cout << "segmentation analysis completed" << endl;
}

void SegmentationManager::load_contig_table()
{
    cout << "loading contig table from " << ifn_contig_table << endl;
    
    ifstream file(ifn_contig_table);
    if (!file.is_open()) {
        cerr << "error: cannot open contig table file " << ifn_contig_table << endl;
        exit(1);
    }
    
    string line;
    int contig_col = -1, length_col = -1;
    
    // parse header
    if (getline(file, line)) {
        istringstream header_iss(line);
        string header;
        int col_index = 0;
        while (getline(header_iss, header, '\t')) {
            if (header == "contig") {
                contig_col = col_index;
            } else if (header == "length") {
                length_col = col_index;
            }
            col_index++;
        }
        
        if (contig_col == -1) {
            cerr << "error: column 'contig' not found in " << ifn_contig_table << endl;
            exit(1);
        }
        if (length_col == -1) {
            cerr << "error: column 'length' not found in " << ifn_contig_table << endl;
            exit(1);
        }
    }
    
    // parse data rows
    contig_lengths.clear();
    while (getline(file, line)) {
        if (line.empty()) continue;
        
        istringstream iss(line);
        string field;
        vector<string> fields;
        while (getline(iss, field, '\t')) {
            fields.push_back(field);
        }
        
        if (static_cast<int>(fields.size()) <= max(contig_col, length_col)) {
            cerr << "error: insufficient columns in line: " << line << endl;
            exit(1);
        }
        
        string contig_id = fields[contig_col];
        uint32_t length = static_cast<uint32_t>(stoul(fields[length_col]));
        contig_lengths[contig_id] = length;
    }
    
    file.close();
    cout << "loaded " << contig_lengths.size() << " contigs from table" << endl;
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
            agg_bp.lib_support[bp.lib_id]++;
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

void SegmentationManager::add_artificial_coords(vector<uint32_t>& coords, uint32_t contig_length)
{
    // skip splitting if max_segment_length is 0
    if (max_segment_length == 0) {
        return;
    }
    
    vector<uint32_t> augmented_coords;
    
    // handle gap from contig start to first coordinate
    uint32_t prev_coord = 1;
    
    for (uint32_t coord : coords) {
        uint32_t gap_length = coord - prev_coord;
        
        // check if gap requires splitting
        if (gap_length > static_cast<uint32_t>(max_segment_length)) {
            // calculate number of segments needed
            int num_segments = static_cast<int>((gap_length + max_segment_length - 1) / max_segment_length);
            
            // calculate even segment size
            double segment_size = static_cast<double>(gap_length) / num_segments;
            
            // add artificial breakpoints
            for (int i = 1; i < num_segments; ++i) {
                uint32_t artificial_coord = prev_coord + static_cast<uint32_t>(i * segment_size);
                augmented_coords.push_back(artificial_coord);
            }
        }
        
        augmented_coords.push_back(coord);
        prev_coord = coord;
    }
    
    // handle gap from last coordinate to contig end
    uint32_t last_coord = coords.empty() ? 1 : coords.back();
    uint32_t final_gap = contig_length - last_coord;
    
    if (final_gap > static_cast<uint32_t>(max_segment_length)) {
        int num_segments = static_cast<int>((final_gap + max_segment_length - 1) / max_segment_length);
        double segment_size = static_cast<double>(final_gap) / num_segments;
        
        for (int i = 1; i < num_segments; ++i) {
            uint32_t artificial_coord = last_coord + static_cast<uint32_t>(i * segment_size);
            augmented_coords.push_back(artificial_coord);
        }
    }
    
    coords = augmented_coords;
}

void SegmentationManager::create_segments_from_coords(const vector<uint32_t>& coords,
                                                      const string& contig_id,
                                                      uint32_t contig_length,
                                                      int& segment_counter)
{
    if (coords.empty()) {
        // no breakpoints, create single segment for entire contig
        string seg_id = "s" + to_string(segment_counter++);
        Segment seg(seg_id, contig_id, 1, contig_length);
        segments.push_back(seg);
        return;
    }
    
    // create segment from start to first coordinate
    if (coords[0] > 1) {
        string seg_id = "s" + to_string(segment_counter++);
        Segment seg(seg_id, contig_id, 1, coords[0]);
        segments.push_back(seg);
    }
    
    // create segments between consecutive coordinates
    for (size_t i = 0; i + 1 < coords.size(); ++i) {
        if (coords[i] == coords[i + 1]) {
            continue;
        }
        string seg_id = "s" + to_string(segment_counter++);
        uint32_t seg_start = coords[i] + 1;
        uint32_t seg_end = coords[i + 1];
        Segment seg(seg_id, contig_id, seg_start, seg_end);
        segments.push_back(seg);
    }
    
    // create segment from last coordinate to end
    if (coords.back() < contig_length) {
        string seg_id = "s" + to_string(segment_counter++);
        uint32_t seg_start = coords.back() + 1;
        Segment seg(seg_id, contig_id, seg_start, contig_length);
        segments.push_back(seg);
    }
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
    
    // group selected breakpoints by contig
    map<string, vector<AggregateBreakpoint*>> contig_breakpoints;
    for (AggregateBreakpoint* bp : selected_bps) {
        contig_breakpoints[bp->contig_id].push_back(bp);
    }
    
    int segment_counter = 1;
    
    // process each contig from loaded contig table
    for (const auto& entry : contig_lengths) {
        const string& contig_id = entry.first;
        uint32_t contig_length = entry.second;
        
        // check if this contig has selected breakpoints
        auto bp_it = contig_breakpoints.find(contig_id);
        
        if (bp_it == contig_breakpoints.end()) {
            // no selected breakpoints for this contig - create single segment
            string seg_id = "s" + to_string(segment_counter++);
            Segment seg(seg_id, contig_id, 1, contig_length);
            segments.push_back(seg);
            continue;
        }
        
        // contig has selected breakpoints - process them
        vector<AggregateBreakpoint*>& bps = bp_it->second;
        
        // sort by coordinate
        sort(bps.begin(), bps.end(),
             [](const AggregateBreakpoint* a, const AggregateBreakpoint* b) {
                 return a->coord < b->coord;
             });
        
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
        
        // extract coordinates from filtered breakpoints
        vector<uint32_t> coords;
        coords.reserve(filtered_bps.size());
        for (const AggregateBreakpoint* bp : filtered_bps) {
            coords.push_back(bp->coord);
        }
        
        // add artificial coordinates where gaps exceed max_segment_length
        add_artificial_coords(coords, contig_length);
        
        // create segments from all coordinates
        create_segments_from_coords(coords, contig_id, contig_length, segment_counter);
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
        const vector<ReadBreakpoint>& breakpoints = entry.second;
        
        for (const ReadBreakpoint& bp : breakpoints) {
            file << bp.lib_id << "\t"
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

