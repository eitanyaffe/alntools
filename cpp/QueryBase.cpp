#include "QueryBase.h"
#include <iostream>

using namespace std;

QueryBase::QueryBase(
    const std::vector<Interval>& intervals,
    ClipMode clip_mode,
    int clip_margin,
    double min_mutations_percent,
    double max_mutations_percent,
    int min_alignment_length,
    int max_alignment_length)
    : intervals(intervals)
    , clip_mode(clip_mode)
    , clip_margin(clip_margin)
    , min_mutations_percent(min_mutations_percent)
    , max_mutations_percent(max_mutations_percent)
    , min_alignment_length(min_alignment_length)
    , max_alignment_length(max_alignment_length)
{
  // constructor implementation - basic initialization done via initializer list
}

bool QueryBase::passes_filter(const Alignment& aln, const AlignmentStore& store) const
{
  return passes_alignment_filter(aln, store, clip_mode, clip_margin, 
                                min_mutations_percent, max_mutations_percent, 
                                min_alignment_length, max_alignment_length);
}
