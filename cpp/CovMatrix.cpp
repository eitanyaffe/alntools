#include "CovMatrix.h"
#include "utils.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>

using namespace std;

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
            
        SegmentInfo seg;
        seg.id = fields[seg_ind];
        seg.contig = fields[contig_ind];
        seg.start = atoi(fields[start_ind].c_str());
        seg.end = atoi(fields[end_ind].c_str());
        seg.length = seg.end - seg.start + 1;
        
        segments.push_back(seg);
    }
    
    cout << "loaded " << segments.size() << " segments" << endl;
}

void CovMatrix::load_libraries(const string& ifn_libraries)
{
    cout << "reading library table: " << ifn_libraries << endl;
    library_files = read_library_table(ifn_libraries);
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

void CovMatrix::calculate_coverage(const SegmentInfo& segment, 
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
    
    // count unique reads
    set<uint32_t> unique_reads;
    for (const auto& alignment_ref : alignments) {
        const Alignment& aln = alignment_ref.get();
        unique_reads.insert(aln.read_index);
    }
    
    int read_count = static_cast<int>(unique_reads.size());
    double length = static_cast<double>(segment.length);
    coverage = read_count / length;
    variance = read_count / (length * length);
}

void CovMatrix::write_fasta(const string& ofn_fasta) const
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
        
        massert(contig_sequences.find(seg.contig) != contig_sequences.end(),
               "contig %s not found in fasta file", seg.contig.c_str());
        
        const string& contig_seq = contig_sequences.at(seg.contig);
        
        // extract subsequence (1-based coordinates, inclusive)
        uint32_t start_0based = seg.start - 1;
        uint32_t length = seg.length;
        
        massert(start_0based + length <= contig_seq.length(),
               "segment %s coordinates out of range for contig %s",
               seg.id.c_str(), seg.contig.c_str());
        
        string seg_seq = contig_seq.substr(start_0based, length);
        
        out << ">" << seg.id << endl;
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
        row.resize(library_files.size(), {0.0, 0.0});
    }
    
    // load each library once and process all segments
    size_t lib_idx = 0;
    for (const auto& lib_entry : library_files) {
        const string& lib_id = lib_entry.first;
        const string& aln_fn = lib_entry.second;
        
        cout << "loading library " << lib_id << " (" << (lib_idx+1) << "/" << library_files.size() << ")" << endl;
        AlignmentStore store;
        store.load(aln_fn);
        
        for (size_t seg_idx = 0; seg_idx < segments.size(); ++seg_idx) {
            const auto& seg = segments[seg_idx];
            
            double cov, var;
            calculate_coverage(seg, store, cov, var);
            
            coverage_matrix[seg_idx][lib_idx] = {cov, var};
            
            if ((seg_idx + 1) % 1000 == 0) {
                cout << "processed " << (seg_idx+1) << "/" << segments.size() << " segments" << endl;
            }
        }
        
        lib_idx++;
    }
    
    // write output file
    ofstream out(ofn_mat.c_str());
    massert(out.is_open(), "could not open file %s", ofn_mat.c_str());
    
    // set output precision
    out << fixed << setprecision(3);
    
    // header
    out << "segment";
    for (size_t i = 0; i < library_files.size(); ++i) {
        out << "\tcov_" << (i+1) << "\tvar_" << (i+1);
    }
    out << endl;
    
    // write data rows
    for (size_t seg_idx = 0; seg_idx < segments.size(); ++seg_idx) {
        out << segments[seg_idx].id;
        for (size_t lib_idx = 0; lib_idx < library_files.size(); ++lib_idx) {
            out << "\t" << coverage_matrix[seg_idx][lib_idx].first 
                << "\t" << coverage_matrix[seg_idx][lib_idx].second;
        }
        out << endl;
    }
    
    out.close();
}

void CovMatrix::compute(const string& ifn_libraries,
                       const string& ifn_segments,
                       const string& ifn_fasta,
                       const string& ofn_mat,
                       const string& ofn_fasta)
{
    load_segments(ifn_segments);
    load_libraries(ifn_libraries);
    load_fasta(ifn_fasta);
    
    write_fasta(ofn_fasta);
    write_matrix(ofn_mat);
    
    cout << "coverage matrix computation complete" << endl;
}

