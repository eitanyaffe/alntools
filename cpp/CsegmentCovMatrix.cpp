#include "CsegmentCovMatrix.h"
#include "utils.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>

using namespace std;

CsegmentCovMatrix::CsegmentCovMatrix() {}

void CsegmentCovMatrix::load_segments(const string& ifn_segments)
{
    cout << "reading segment table: " << ifn_segments << endl;
    ifstream file(ifn_segments);
    massert(file.is_open(), "could not open file %s", ifn_segments.c_str());
    
    string line;
    int seg_ind = -1, contig_ind = -1, start_ind = -1, end_ind = -1;
    
    // parse header
    if (getline(file, line)) {
        istringstream header_iss(line);
        string header;
        int col_index = 0;
        while (getline(header_iss, header, '\t')) {
            if (header == "segment") seg_ind = col_index;
            else if (header == "contig") contig_ind = col_index;
            else if (header == "start") start_ind = col_index;
            else if (header == "end") end_ind = col_index;
            col_index++;
        }
        
        massert(seg_ind >= 0, "column 'segment' not found in %s", ifn_segments.c_str());
        massert(contig_ind >= 0, "column 'contig' not found in %s", ifn_segments.c_str());
        massert(start_ind >= 0, "column 'start' not found in %s", ifn_segments.c_str());
        massert(end_ind >= 0, "column 'end' not found in %s", ifn_segments.c_str());
    }
    
    int total_count = 0;
    int filtered_count = 0;
    
    // parse data rows
    while (getline(file, line)) {
        if (line.empty()) continue;
        
        istringstream iss(line);
        vector<string> fields;
        string field;
        while (getline(iss, field, '\t')) {
            fields.push_back(field);
        }
        
        if (fields.size() == 0) continue;
            
        CsegSegment seg;
        seg.id = fields[seg_ind];
        seg.contig = fields[contig_ind];
        seg.start = atoi(fields[start_ind].c_str());
        seg.end = atoi(fields[end_ind].c_str());
        seg.length = seg.end - seg.start + 1;
        seg.index = segments.size();
        
        total_count++;
        
        // filter by minimum segment length
        if (seg.length < min_segment_length) {
            filtered_count++;
            continue;
        }
        
        segments.push_back(seg);
    }
    
    cout << "loaded " << segments.size() << " segments" << endl;
    if (filtered_count > 0) {
        cout << "filtered " << filtered_count << " segments shorter than " << min_segment_length << " bp" << endl;
    }
}

void CsegmentCovMatrix::load_libraries(const string& ifn_libraries)
{
    cout << "reading library table: " << ifn_libraries << endl;
    ifstream file(ifn_libraries);
    massert(file.is_open(), "could not open file %s", ifn_libraries.c_str());
    
    string line;
    int lib_id_ind = -1, aln_fn_ind = -1;
    
    // parse header
    if (getline(file, line)) {
        istringstream header_iss(line);
        string header;
        int col_index = 0;
        while (getline(header_iss, header, '\t')) {
            if (header == "lib_id") lib_id_ind = col_index;
            else if (header == "aln_fn") aln_fn_ind = col_index;
            col_index++;
        }
        
        massert(lib_id_ind >= 0, "column 'lib_id' not found in %s", ifn_libraries.c_str());
        massert(aln_fn_ind >= 0, "column 'aln_fn' not found in %s", ifn_libraries.c_str());
    }
    
    // parse data rows
    while (getline(file, line)) {
        if (line.empty()) continue;
        
        istringstream iss(line);
        vector<string> fields;
        string field;
        while (getline(iss, field, '\t')) {
            fields.push_back(field);
        }
        
        if (fields.size() == 0) continue;
        
        string lib_id = fields[lib_id_ind];
        string aln_fn = fields[aln_fn_ind];
        
        library_files[lib_id] = aln_fn;
        library_ids.push_back(lib_id);
    }
    
    cout << "loaded " << library_ids.size() << " libraries" << endl;
}

void CsegmentCovMatrix::load_cseg_map(const string& ifn_cseg_map)
{
    cout << "reading cseg map: " << ifn_cseg_map << endl;
    ifstream file(ifn_cseg_map);
    massert(file.is_open(), "could not open file %s", ifn_cseg_map.c_str());
    
    string line;
    int segment_ind = -1, csegment_ind = -1;
    
    // parse header
    if (getline(file, line)) {
        istringstream header_iss(line);
        string header;
        int col_index = 0;
        while (getline(header_iss, header, '\t')) {
            if (header == "segment") segment_ind = col_index;
            else if (header == "csegment") csegment_ind = col_index;
            col_index++;
        }
        
        massert(segment_ind >= 0, "column 'segment' not found in %s", ifn_cseg_map.c_str());
        massert(csegment_ind >= 0, "column 'csegment' not found in %s", ifn_cseg_map.c_str());
    }
    
    // build segment ID to index map
    unordered_map<string, size_t> segment_id_to_index;
    for (size_t i = 0; i < segments.size(); ++i) {
        segment_id_to_index[segments[i].id] = i;
    }
    
    // parse data rows
    while (getline(file, line)) {
        if (line.empty()) continue;
        
        istringstream iss(line);
        vector<string> fields;
        string field;
        while (getline(iss, field, '\t')) {
            fields.push_back(field);
        }
        
        if (fields.size() == 0) continue;
        
        string segment_id = fields[segment_ind];
        string csegment_id = fields[csegment_ind];
        
        auto it = segment_id_to_index.find(segment_id);
        if (it == segment_id_to_index.end()) {
            continue;
        }
        
        size_t seg_idx = it->second;
        csegment_to_segment_indices[csegment_id].push_back(seg_idx);
        cseg_to_cluster[segment_id] = csegment_id;
    }
    
    cout << "loaded " << csegment_to_segment_indices.size() << " csegments" << endl;
}

void CsegmentCovMatrix::count_total_bp_per_csegment(const string& csegment,
                                                     const AlignmentStore& store,
                                                     uint64_t& total_bp) const
{
    total_bp = 0;
    
    auto it = csegment_to_segment_indices.find(csegment);
    if (it == csegment_to_segment_indices.end()) {
        return;
    }
    
    const vector<size_t>& segment_indices = it->second;
    
    // map from read_index to maximum intersection length for that read
    unordered_map<uint32_t, uint32_t> read_to_max_bp;
    
    for (size_t seg_idx : segment_indices) {
        const CsegSegment& seg = segments[seg_idx];
        
        if (!store.has_contig_id(seg.contig)) {
            continue;
        }
        
        // convert 1-based inclusive to 0-based half-open
        uint32_t start_0based = seg.start - 1;
        uint32_t end_0based = seg.end;
        
        Interval interval(seg.contig, start_0based, end_0based);
        auto alignments = store.get_alignments_intersecting_interval(interval);
        
        for (const auto& alignment_ref : alignments) {
            const Alignment& aln = alignment_ref.get();
            
            // apply alignment filter
            if (!passes_alignment_filter(aln, store, clip_mode, clip_margin,
                                        min_mutations_percent, max_mutations_percent,
                                        min_alignment_length, max_alignment_length)) {
                continue;
            }
            
            // calculate intersection length
            uint32_t intersection_start = max(start_0based, aln.contig_start);
            uint32_t intersection_end = min(end_0based, aln.contig_end);
            
            if (intersection_end > intersection_start) {
                uint32_t intersection_bp = intersection_end - intersection_start;
                
                // update map to keep maximum intersection length per read
                auto read_it = read_to_max_bp.find(aln.read_index);
                if (read_it == read_to_max_bp.end()) {
                    read_to_max_bp[aln.read_index] = intersection_bp;
                } else {
                    if (intersection_bp > read_it->second) {
                        read_it->second = intersection_bp;
                    }
                }
            }
        }
    }
    
    // sum all maximum intersection lengths (each read contributes once with its max)
    for (const auto& entry : read_to_max_bp) {
        total_bp += entry.second;
    }
}

void CsegmentCovMatrix::write_coverage(const string& ofn_coverage)
{
    cout << "writing coverage matrix: " << ofn_coverage << endl;
    
    // store coverage matrix: csegment -> library -> total_bp
    map<string, map<string, uint64_t>> coverage_matrix;
    
    // initialize all csegments
    for (const auto& entry : csegment_to_segment_indices) {
        coverage_matrix[entry.first] = map<string, uint64_t>();
    }
    
    // process each library
    for (size_t lib_idx = 0; lib_idx < library_ids.size(); ++lib_idx) {
        const string& lib_id = library_ids[lib_idx];
        auto it = library_files.find(lib_id);
        if (it == library_files.end()) {
            continue;
        }
        const string& aln_fn = it->second;
        
        cout << "loading library " << lib_id << " (" << (lib_idx+1) << "/" << library_ids.size() << ")" << endl;
        AlignmentStore store;
        store.load(aln_fn);
        cout << "loaded " << store.get_read_count() << " reads, " << store.get_alignment_count() << " alignments" << endl;
        
        // initialize short indel counting for mutation density calculations
        store.count_short_indels(min_indel_length);
        
        // initialize read alignment index if LOCAL_ALIGN mode is used
        init_local_align_if_needed(store, clip_mode);
        cout << "processing csegments..." << endl;
        
        // process each csegment
        size_t total_csegs = csegment_to_segment_indices.size();
        size_t cseg_idx = 0;
        for (const auto& entry : csegment_to_segment_indices) {
            const string& csegment = entry.first;
            uint64_t total_bp = 0;
            
            count_total_bp_per_csegment(csegment, store, total_bp);
            
            coverage_matrix[csegment][lib_id] = total_bp;
            
            cseg_idx++;
            if (cseg_idx % 1000 == 0) {
                cout << "processed " << cseg_idx << "/" << total_csegs << " csegments for library " << lib_id << endl;
            }
        }
        cout << "completed processing " << total_csegs << " csegments for library " << lib_id << endl;
    }
    
    // write output file
    ofstream out(ofn_coverage.c_str());
    massert(out.is_open(), "could not open file %s", ofn_coverage.c_str());
    
    // header
    out << "csegment";
    for (const string& lib_id : library_ids) {
        out << "\t" << lib_id;
    }
    out << endl;
    
    // write data rows
    for (const auto& cseg_entry : coverage_matrix) {
        out << cseg_entry.first;
        for (const string& lib_id : library_ids) {
            auto lib_it = cseg_entry.second.find(lib_id);
            if (lib_it != cseg_entry.second.end()) {
                out << "\t" << lib_it->second;
            } else {
                out << "\t0";
            }
        }
        out << endl;
    }
    
    out.close();
}

void CsegmentCovMatrix::write_read_counts(const string& ofn_read_counts,
                                          const map<string, uint32_t>& lib_total_reads) const
{
    cout << "writing read counts: " << ofn_read_counts << endl;
    ofstream out(ofn_read_counts.c_str());
    massert(out.is_open(), "could not open file %s", ofn_read_counts.c_str());
    
    out << "lib_id\ttotal_reads" << endl;
    for (const auto& entry : lib_total_reads) {
        out << entry.first << "\t" << entry.second << endl;
    }
    
    out.close();
}

void CsegmentCovMatrix::compute(const string& ifn_libraries,
                               const string& ifn_segments,
                               const string& ifn_cseg_map,
                               const string& ofn_coverage,
                               const string& ofn_read_counts,
                               uint32_t min_seg_len,
                               ClipMode clip_mode_param,
                               int clip_margin_param,
                               double min_mutations_percent_param,
                               double max_mutations_percent_param,
                               int min_alignment_length_param,
                               int max_alignment_length_param,
                               int min_indel_length_param)
{
    min_segment_length = min_seg_len;
    clip_mode = clip_mode_param;
    clip_margin = clip_margin_param;
    min_mutations_percent = min_mutations_percent_param;
    max_mutations_percent = max_mutations_percent_param;
    min_alignment_length = min_alignment_length_param;
    max_alignment_length = max_alignment_length_param;
    min_indel_length = min_indel_length_param;
    
    load_segments(ifn_segments);
    load_libraries(ifn_libraries);
    load_cseg_map(ifn_cseg_map);
    
    cout << "collecting total read counts per library..." << endl;
    // collect total read counts per library
    map<string, uint32_t> lib_total_reads;
    for (const string& lib_id : library_ids) {
        auto it = library_files.find(lib_id);
        if (it == library_files.end()) {
            continue;
        }
        AlignmentStore store;
        store.load(it->second);
        lib_total_reads[lib_id] = store.get_read_count();
        cout << "library " << lib_id << ": " << lib_total_reads[lib_id] << " total reads" << endl;
    }
    
    cout << "processing coverage for " << csegment_to_segment_indices.size() << " csegments across " << library_ids.size() << " libraries..." << endl;
    write_coverage(ofn_coverage);
    write_read_counts(ofn_read_counts, lib_total_reads);
    
    cout << "csegment coverage computation complete" << endl;
}

