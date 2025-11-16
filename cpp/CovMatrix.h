#pragma once

#include "alignment_store.h"
#include "utils.h"
#include <string>
#include <vector>
#include <map>
#include <unordered_map>

struct CovSegment {
    std::string id;
    std::string contig;
    uint32_t start;  // 1-based inclusive
    uint32_t end;    // 1-based inclusive
    uint32_t length;
    size_t index;

    CovSegment() : start(0), end(0), length(0), index(0) {}
};

class CovMatrix {
private:
    std::vector<CovSegment> segments;
    map<string, string> library_files;
    vector<string> library_ids;
    unordered_map<string, string> contig_sequences;
    uint32_t min_segment_length;
    
    ClipMode clip_mode;
    int clip_margin;
    double min_mutations_percent;
    double max_mutations_percent;
    int min_alignment_length;
    int max_alignment_length;
    int min_indel_length;
    
    void load_segments(const string& ifn_segments);
    void load_libraries(const string& ifn_libraries);
    void load_fasta(const std::string& ifn_fasta);
    
    void calculate_coverage(const CovSegment& segment, 
                          const AlignmentStore& store,
                          double& coverage, 
                          double& variance) const;
    
    void write_fasta(const string& ofn_fasta, bool actual_nts) const;
    void write_matrix(const string& ofn_mat);
    void write_lib_map(const string& ofn_lib_map) const;

public:
    CovMatrix();
    
    void compute(const string& ifn_libraries,
                const string& ifn_segments,
                const string& ifn_fasta,
                const string& ofn_mat,
                const string& ofn_fasta,
                bool actual_nts,
                bool should_create_fasta,
                uint32_t min_segment_length,
                ClipMode clip_mode,
                int clip_margin,
                double min_mutations_percent,
                double max_mutations_percent,
                int min_alignment_length,
                int max_alignment_length,
                int min_indel_length,
                const string& ofn_lib_map = "");
};

