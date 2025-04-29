#ifndef QUERYBIN_H
#define QUERYBIN_H

#include "alignment_store.h" // Includes aln_types.h indirectly
#include <map>
#include <string>
#include <utility> // For std::pair
#include <vector>

// Data structure to hold aggregated results for a single bin
struct BinData {
  int sequenced_basepairs = 0;
  int mutation_count = 0;
};

// Data structure representing a single row in the bin output file
struct BinOutputRow {
  std::string contig;
  uint32_t bin_start;
  uint32_t bin_end;
  int bin_length;
  int sequenced_basepairs;
  int mutation_count;
};

class QueryBin {
  private:
  const std::vector<Interval>& intervals;
  const AlignmentStore& store;
  int binsize;

  // Use a map to store results, keyed by {contig_index, bin_start}
  std::map<std::pair<uint32_t, uint32_t>, BinData> bin_results;
  // Vector to store the formatted output rows before writing
  std::vector<BinOutputRow> output_rows;

  // execute the query
  void aggregate_data();

  // generate the output rows
  void generate_output_rows();

  public:
  QueryBin(
      const std::vector<Interval>& intervals,
      const AlignmentStore& store,
      int binsize);

  // execute the query
  void execute();

  // write the output rows to a table
  void write_to_csv(const std::string& ofn_prefix);

  // Getter for R interface
  const std::vector<BinOutputRow>& get_output_rows() const { return output_rows; }
};

#endif // QUERYBIN_H