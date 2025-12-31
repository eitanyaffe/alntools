#pragma once

#include "alignment_store.h"
#include "utils.h"
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>

struct CsegSegment {
    std::string id;
    std::string contig;
    uint32_t start;  // 1-based inclusive
    uint32_t end;    // 1-based inclusive
    uint32_t length;
    size_t index;

    CsegSegment() : start(0), end(0), length(0), index(0) {}
};

class CsegmentCovMatrix {
private:
    std::vector<CsegSegment> segments;
    std::map<std::string, std::string> library_files;
    std::vector<std::string> library_ids;
    std::unordered_map<std::string, std::string> cseg_to_cluster;
    std::unordered_map<std::string, std::vector<size_t>> csegment_to_segment_indices;
    
    uint32_t min_segment_length;
    ClipMode clip_mode;
    int clip_margin;
    double min_mutations_percent;
    double max_mutations_percent;
    int min_alignment_length;
    int max_alignment_length;
    int min_indel_length;
    
    void load_segments(const std::string& ifn_segments);
    void load_query_table(const std::string& ifn_query);
    void load_libraries(const std::string& ifn_libraries);
    void load_cseg_map(const std::string& ifn_cseg_map);
    
    void count_unique_reads_per_csegment(const std::string& csegment,
                                        const AlignmentStore& store,
                                        uint32_t& unique_read_count) const;
    
    void write_coverage(const std::string& ofn_coverage);
    void write_read_counts(const std::string& ofn_read_counts,
                          const std::map<std::string, uint32_t>& lib_total_reads) const;

public:
    CsegmentCovMatrix();
    
    void compute(const std::string& ifn_libraries,
                const std::string& ifn_segments,
                const std::string& ifn_cseg_map,
                const std::string& ofn_coverage,
                const std::string& ofn_read_counts,
                uint32_t min_segment_length,
                ClipMode clip_mode,
                int clip_margin,
                double min_mutations_percent,
                double max_mutations_percent,
                int min_alignment_length,
                int max_alignment_length,
                int min_indel_length);
};

