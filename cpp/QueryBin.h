#ifndef QUERYBIN_H
#define QUERYBIN_H

#include "alignment_store.h" // Includes aln_types.h indirectly
#include <map>
#include <string>
#include <utility> // For std::pair
#include <vector>

// Data structure to hold aggregated results for a single bin
struct BinData {
  long long sequenced_basepairs = 0;
  int mutation_count = 0;
};

class QueryBin {
  private:
  const std::vector<Interval>& intervals;
  std::string ofn_prefix;
  const AlignmentStore& store;
  int binsize;

  // Use a map to store results, keyed by {contig_index, bin_start}
  std::map<std::pair<uint32_t, uint32_t>, BinData> bin_results;
  bool skip_empty_bins;
  void aggregate_data(); // Helper function to perform the aggregation
  void write_data_to_file(); // Helper function to write the aggregated data

  public:
  QueryBin(
      const std::vector<Interval>& intervals,
      const std::string& ofn_prefix,
      const AlignmentStore& store,
      int binsize,
      bool skip_empty_bins);
  void write_to_csv(); // Main public method to coordinate aggregation and writing
};

#endif // QUERYBIN_H