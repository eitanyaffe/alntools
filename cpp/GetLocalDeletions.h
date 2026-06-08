#pragma once

#include "alignment_store.h"
#include "utils.h"
#include <fstream>
#include <string>
#include <vector>

class GetLocalDeletions {
private:
    int margin;
    int flank;
    int min_indel_length;

    // format non-short-indel mutations in [range_start, range_end] as "pos[desc],..." or "."
    std::string format_flanking_mutations(const Alignment& aln,
                                          const AlignmentStore& store,
                                          uint32_t range_start,
                                          uint32_t range_end) const;

    void find_deletions(const AlignmentStore& store,
                        const std::vector<Interval>& intervals,
                        std::ofstream& out,
                        bool emit_interval_id) const;

public:
    GetLocalDeletions();

    void compute(const std::string& ifn_aln,
                 const std::string& ifn_intervals,
                 const std::string& ofn,
                 int margin,
                 int flank,
                 int min_indel_length);
};
