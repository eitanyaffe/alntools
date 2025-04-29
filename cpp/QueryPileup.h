#ifndef QUERYPILEUP_H
#define QUERYPILEUP_H

#include "alignment_store.h" // Includes aln_types.h indirectly
#include <map>
#include <string>
#include <utility> // For std::pair
#include <vector>

// Enum to control output verbosity for pileup
enum class PileupReportMode {
  ALL, // Report every position within the query intervals
  COVERED, // Report only positions with coverage > 0
  MUTATED // Report only positions with at least one mutation observed
};

// Data structure to hold aggregated results for a single genomic position
struct PileupData {
  int coverage = 0;
  // Stores counts for specific mutations observed at this position
  std::map<std::string, int> mutation_counts; // Key: Mutation::to_string(), Value: count
};

// Data structure representing a single row in the pileup output file
struct PileupOutputRow {
  std::string contig;
  uint32_t position; // 1-based
  std::string variant;
  int count;
  int coverage;
  int cumsum;
};

class QueryPileup {
  private:
  const std::vector<Interval>& intervals;
  const AlignmentStore& store;
  PileupReportMode report_mode;

  // Use a map to store results, keyed by {contig_index, position (0-based)}
  std::map<std::pair<uint32_t, uint32_t>, PileupData> pileup_results;
  // Vector to store the formatted output rows before writing
  std::vector<PileupOutputRow> output_rows;

  void aggregate_data();
  void generate_output_rows();
  void write_rows_to_file(const std::string& ofn_prefix);

  public:
  QueryPileup(const std::vector<Interval>& intervals, const AlignmentStore& store, PileupReportMode report_mode);

  // execute the query
  void execute();

  // write the output rows to a table
  void write_to_csv(const std::string& ofn_prefix);
};

#endif // QUERYPILEUP_H