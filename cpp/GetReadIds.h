#pragma once

#include "alignment_store.h"
#include "utils.h"
#include <string>
#include <unordered_map>
#include <vector>

struct ReadSegment {
    std::string contig;
    uint32_t start;  // 1-based inclusive
    uint32_t end;    // 1-based inclusive
    std::string bin_id;
};

// best bin assignment for a read: the bin whose segment has the longest intersection
struct ReadBinEntry {
    std::string bin_id;
    uint32_t length;

    ReadBinEntry() : length(0) {}
    ReadBinEntry(const std::string& bin_id, uint32_t length) : bin_id(bin_id), length(length) {}
};

class GetReadIds {
private:
    std::vector<ReadSegment> segments;

    ClipMode clip_mode;
    int clip_margin;
    double min_mutations_percent;
    double max_mutations_percent;
    int min_alignment_length;
    int max_alignment_length;
    int min_indel_length;

    void load_segments(const std::string& ifn_segments);
    void collect_read_ids(AlignmentStore& store,
                          std::unordered_map<uint32_t, ReadBinEntry>& read_map) const;
    void write_output(const std::string& ofn,
                      const AlignmentStore& store,
                      const std::unordered_map<uint32_t, ReadBinEntry>& read_map) const;

public:
    GetReadIds();

    void compute(const std::string& ifn_aln,
                 const std::string& ifn_segments,
                 const std::string& ofn,
                 ClipMode clip_mode,
                 int clip_margin,
                 double min_mutations_percent,
                 double max_mutations_percent,
                 int min_alignment_length,
                 int max_alignment_length,
                 int min_indel_length);
};
