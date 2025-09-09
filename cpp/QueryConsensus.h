#ifndef QUERYCONSENSUS_H
#define QUERYCONSENSUS_H

#include "QueryBase.h"
#include <map>
#include <set>
#include <string>
#include <utility> // For std::pair
#include <vector>

// Data structure to hold variant information for a single position
struct VariantData {
  int count = 0;
  int coverage = 0;
  
  VariantData() = default;
  VariantData(int c, int cov) : count(c), coverage(cov) {}
};

// Data structure representing a single row in the consensus output file
struct ConsensusOutputRow {
  std::string contig;
  uint32_t position;
  std::string variant_type; // SUB, INS, DEL
  std::string variant_desc; // nucleotide description
  int count;
  int coverage;
  double frequency;
};

class QueryConsensus : public QueryBase {
  private:
  const AlignmentStore& store;
  double consensus_threshold;
  int min_consensus_coverage;
  int num_threads;

  // Use a map to store results, keyed by {contig_index, position, variant_key}
  std::map<std::tuple<uint32_t, uint32_t, std::string>, VariantData> variant_results;
  // Vector to store the formatted output rows before writing
  std::vector<ConsensusOutputRow> output_rows;

  // build contig to intervals mapping for efficient lookup
  void build_contig_to_intervals_map();

  // process a single alignment into variant results (thread-safe)
  void process_single_alignment(const Alignment& aln, std::map<std::tuple<uint32_t, uint32_t, std::string>, VariantData>& target_variant_results);

  // merge thread-local VariantData into global variant_results
  void merge_variant_data(const std::map<std::tuple<uint32_t, uint32_t, std::string>, VariantData>& local_data);

  // calculate coverage for all variant positions
  void calc_variant_coverage();

  // generate the output rows from variant results
  void generate_output_rows();

  public:
  QueryConsensus(
      const std::vector<Interval>& intervals,
      const AlignmentStore& store,
      double consensus_threshold = 0.9,
      int min_consensus_coverage = 5,
      int num_threads = 0,
      ClipMode clip_mode = ClipMode::ALL,
      int clip_margin = 10,
      double min_mutations_percent = 0.0,
      double max_mutations_percent = 10.0,
      int min_alignment_length = 0,
      int max_alignment_length = 0);

  // execute the query
  void execute();
  
  // execute data aggregation
  void aggregate_data();

  // write the output rows to a table
  void write_to_csv(const std::string& ofn_prefix);

  // Getter for R interface
  const std::vector<ConsensusOutputRow>& get_output_rows() const { return output_rows; }
};

#endif // QUERYCONSENSUS_H
