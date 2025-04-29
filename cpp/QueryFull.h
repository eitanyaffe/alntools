#ifndef QUERYFULL_H
#define QUERYFULL_H

#include "alignment_store.h"
#include <string>
#include <vector>

class QueryFull {
  private:
  std::vector<Interval> intervals;
  std::string ofn_prefix;
  const AlignmentStore& store;

  public:
  QueryFull(const std::vector<Interval>& intervals, const std::string& ofn_prefix, const AlignmentStore& store);

  // Write the full query results to a CSV file
  void write_to_csv();
};

#endif // QUERYFULL_H