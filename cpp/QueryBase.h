#ifndef QUERYBASE_H
#define QUERYBASE_H

#include "alignment_store.h" // Includes aln_types.h indirectly
#include "utils.h"
#include <map>
#include <string>
#include <vector>

class QueryBase {
protected:
  const std::vector<Interval>& intervals;
  ClipMode clip_mode;
  int clip_margin;
  double min_mutations_percent;
  double max_mutations_percent;
  int min_alignment_length;
  int max_alignment_length;

  // shared efficiency optimization - each derived class builds this with their own store(s)
  std::map<uint32_t, std::vector<const Interval*>> contig_to_intervals_map;
  
  // store-agnostic utility methods
  bool passes_filter(const Alignment& aln, const AlignmentStore& store) const;

public:
  QueryBase(
      const std::vector<Interval>& intervals,
      ClipMode clip_mode = ClipMode::ALL,
      int clip_margin = 10,
      double min_mutations_percent = 0.0,
      double max_mutations_percent = 10.0,
      int min_alignment_length = 0,
      int max_alignment_length = 0);

  virtual ~QueryBase() = default;
  
  // pure virtual methods that derived classes must implement
  virtual void execute() = 0;
  virtual void write_to_csv(const std::string& ofn_prefix) = 0;
};

#endif // QUERYBASE_H
