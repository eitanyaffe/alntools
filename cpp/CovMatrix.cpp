#include "CovMatrix.h"
#include "utils.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <random>

using namespace std;

namespace {
string generate_random_nts(uint32_t length)
{
    static thread_local std::mt19937 rng(std::random_device{}());
    static const char bases[] = {'A', 'C', 'G', 'T'};
    std::uniform_int_distribution<int> dist(0, 3);

    string nts;
    nts.reserve(length);
    for (uint32_t i = 0; i < length; ++i) {
        nts.push_back(bases[dist(rng)]);
    }
    return nts;
}
}

CovMatrix::CovMatrix() {}

void CovMatrix::load_segments(const string& ifn_segments)
{
    cout << "reading segment table: " << ifn_segments << endl;
    ifstream file(ifn_segments);
    massert(file.is_open(), "could not open file %s", ifn_segments.c_str());
    
    string line;
    vector<string> headers;
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
            
        CovSegment seg;
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

void CovMatrix::load_libraries(const string& ifn_libraries)
{
    cout << "reading library table: " << ifn_libraries << endl;
    read_library_table(ifn_libraries, library_files, library_ids);
}

void CovMatrix::load_fasta(const string& ifn_fasta)
{
    unordered_set<string> all_contigs;
    for (const auto& seg : segments) {
        all_contigs.insert(seg.contig);
    }
    read_fasta(ifn_fasta, all_contigs, contig_sequences);
    cout << "loaded " << contig_sequences.size() << " contigs" << endl;
}

void CovMatrix::calculate_coverage(const CovSegment& segment, 
                                   const AlignmentStore& store,
                                   double& coverage, 
                                   double& variance) const
{
    coverage = 0.0;
    variance = 0.0;
    
    if (!store.has_contig_id(segment.contig)) {
        return;
    }
    
    // convert 1-based inclusive to 0-based half-open
    uint32_t start_0based = segment.start - 1;
    uint32_t end_0based = segment.end;
    
    Interval interval(segment.contig, start_0based, end_0based);
    auto alignments = store.get_alignments_intersecting_interval(interval);
    
    // sum intersection lengths with filtering
    uint32_t total_bp = 0;
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
            total_bp += (intersection_end - intersection_start);
        }
    }
    
    double length = static_cast<double>(segment.length);
    coverage = total_bp / length;
    variance = total_bp / (length * length);
}

void CovMatrix::write_fasta(const string& ofn_fasta, bool actual_nts) const
{
    cout << "writing segment fasta: " << ofn_fasta << endl;
    ofstream out(ofn_fasta.c_str());
    massert(out.is_open(), "could not open file %s", ofn_fasta.c_str());
    
    for (const auto& seg : segments) {
        // validate segment coordinates
        if (seg.start == 0) {
            cerr << "error: segment " << seg.id << " has invalid start coordinate 0 (coordinates should be 1-based)" << endl;
            cerr << "  segment: " << seg.id << ", contig: " << seg.contig 
                 << ", start: " << seg.start << ", end: " << seg.end << ", length: " << seg.length << endl;
            exit(1);
        }
        if (seg.length == 0) {
            cerr << "error: segment " << seg.id << " has zero length" << endl;
            cerr << "  segment: " << seg.id << ", contig: " << seg.contig 
                 << ", start: " << seg.start << ", end: " << seg.end << ", length: " << seg.length << endl;
            exit(1);
        }
        if (seg.start > seg.end) {
            cerr << "error: segment " << seg.id << " has start > end" << endl;
            cerr << "  segment: " << seg.id << ", contig: " << seg.contig 
                 << ", start: " << seg.start << ", end: " << seg.end << ", length: " << seg.length << endl;
            exit(1);
        }
        
        out << ">" << seg.id << endl;

        if (!actual_nts) {
            out << generate_random_nts(seg.length) << endl;
            continue;
        }

        auto contig_it = contig_sequences.find(seg.contig);
        massert(contig_it != contig_sequences.end(),
               "contig %s not found in fasta file", seg.contig.c_str());

        const string& contig_seq = contig_it->second;

        // extract subsequence (1-based coordinates, inclusive)
        uint32_t start_0based = seg.start - 1;
        uint32_t length = seg.length;

        massert(start_0based + length <= contig_seq.length(),
               "segment %s coordinates out of range for contig %s",
               seg.id.c_str(), seg.contig.c_str());

        string seg_seq = contig_seq.substr(start_0based, length);

        out << seg_seq << endl;
    }
    
    out.close();
}

void CovMatrix::write_matrix(const string& ofn_mat)
{
    cout << "writing coverage matrix: " << ofn_mat << endl;
    
    // store coverage matrix in memory: segment_index -> library_index -> (cov, var)
    vector<vector<pair<double, double>>> coverage_matrix(segments.size());
    for (auto& row : coverage_matrix) {
        row.resize(library_ids.size(), {0.0, 0.0});
    }
    
    // load each library once and process all segments in library_ids order
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
        
        // initialize short indel counting for mutation density calculations
        store.count_short_indels(min_indel_length);
        
        // initialize read alignment index if LOCAL_ALIGN mode is used
        init_local_align_if_needed(store, clip_mode);
        
        for (size_t seg_idx = 0; seg_idx < segments.size(); ++seg_idx) {
            const auto& seg = segments[seg_idx];
            
            double cov, var;
            calculate_coverage(seg, store, cov, var);
            
            coverage_matrix[seg_idx][lib_idx] = {cov, var};
            
            if ((seg_idx + 1) % 1000 == 0) {
                cout << "processed " << (seg_idx+1) << "/" << segments.size() << " segments" << endl;
            }
        }
    }
    
    // write output file
    ofstream out(ofn_mat.c_str());
    massert(out.is_open(), "could not open file %s", ofn_mat.c_str());
    
    // set output precision
    out << fixed << setprecision(3);
    
    // header
    out << "segment";
    for (size_t i = 0; i < library_ids.size(); ++i) {
        out << "\tcov_" << (i+1) << "\tvar_" << (i+1);
    }
    out << endl;
    
    // write data rows
    for (size_t seg_idx = 0; seg_idx < segments.size(); ++seg_idx) {
        out << segments[seg_idx].id;
        for (size_t lib_idx = 0; lib_idx < library_ids.size(); ++lib_idx) {
            out << "\t" << coverage_matrix[seg_idx][lib_idx].first 
                << "\t" << coverage_matrix[seg_idx][lib_idx].second;
        }
        out << endl;
    }
    
    out.close();
}

void CovMatrix::write_lib_map(const string& ofn_lib_map) const
{
    cout << "writing library mapping to " << ofn_lib_map << endl;
    ofstream out(ofn_lib_map.c_str());
    massert(out.is_open(), "could not open file %s", ofn_lib_map.c_str());
    
    out << "lib_index\tlib_id" << endl;
    for (size_t i = 0; i < library_ids.size(); ++i) {
        out << "lib_" << (i+1) << "\t" << library_ids[i] << endl;
    }
    
    out.close();
}

void CovMatrix::compute(const string& ifn_libraries,
                       const string& ifn_segments,
                       const string& ifn_fasta,
                       const string& ofn_mat,
                       const string& ofn_fasta,
                       bool actual_nts,
                       bool should_create_fasta,
                       uint32_t min_seg_len,
                       ClipMode clip_mode_param,
                       int clip_margin_param,
                       double min_mutations_percent_param,
                       double max_mutations_percent_param,
                       int min_alignment_length_param,
                       int max_alignment_length_param,
                       int min_indel_length_param,
                       const string& ofn_lib_map)
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

    if (should_create_fasta && actual_nts) {
        load_fasta(ifn_fasta);
    }

    if (should_create_fasta) {
        write_fasta(ofn_fasta, actual_nts);
    }
    
    write_matrix(ofn_mat);
    
    if (!ofn_lib_map.empty()) {
        write_lib_map(ofn_lib_map);
    }
    
    cout << "coverage matrix computation complete" << endl;
}

