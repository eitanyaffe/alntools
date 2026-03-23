#ifndef QUERYBIN_H
#define QUERYBIN_H

#include "QueryBase.h"
#include <map>
#include <set>
#include <string>
#include <utility> // For std::pair
#include <vector>

// Data structure to hold aggregated results for a single bin
struct BinData {
  int sequenced_basepairs = 0;
  int mutation_count = 0;
  std::set<uint32_t> unique_reads; // track unique read indices for read count
  std::map<std::string, int> variant_counts; // key: position_type_nts, value: count
  std::map<std::string, int> position_coverage; // key: position, value: read count at position
  std::map<std::string, int> clip_left_counts; // key: position, value: count of left clips
  std::map<std::string, int> clip_right_counts; // key: position, value: count of right clips
  // mutation distance categories: counts of alignments in each category (per bp)
  int dist_none = 0;         // exactly 0 mutations
  int dist_5 = 0;            // 1e-5 to 1e-4 per bp (10-100 mutations per 100kb)
  int dist_4 = 0;            // 1e-4 to 1e-3 per bp (100-1000 mutations per 100kb)
  int dist_3 = 0;            // 1e-3 to 1e-2 per bp (1000-10000 mutations per 100kb)
  int dist_2 = 0;            // 1e-2 to 1e-1 per bp (10000-100000 mutations per 100kb)
  int dist_1_plus = 0;       // above 1e-1 per bp (>100000 mutations per 100kb)
  // mutation density values for median calculation
  std::vector<double> mutation_densities; // collect local_mutations_per_bp for each alignment
};

// Data structure representing a single row in the bin output file
struct BinOutputRow {
  std::string contig;
  uint32_t bin_start;
  uint32_t bin_end;
  int bin_length;
  int sequenced_basepairs;
  int read_count; // number of unique reads
  int mutation_count;
  double median_mutation_density; // median mutations per bp across alignments in bin
  double seg_sites_density; // segregating sites per bp
  double non_ref_sites_density; // non-reference sites per bp
  double seg_clip_density; // segregating clip sites per bp
  double non_ref_clip_density; // non-reference clip sites per bp
  int dist_none;         // exactly 0 mutations
  int dist_5;            // 1e-5 to 1e-4 per bp
  int dist_4;            // 1e-4 to 1e-3 per bp
  int dist_3;            // 1e-3 to 1e-2 per bp
  int dist_2;            // 1e-2 to 1e-1 per bp
  int dist_1_plus;       // above 1e-1 per bp
};

class QueryBin : public QueryBase {
  private:
  const AlignmentStore& store;
  int binsize;
  double seg_threshold;
  double non_ref_threshold;
  int num_threads;
  int min_seg_support;

  // Use a map to store results, keyed by {contig_index, bin_start}
  std::map<std::pair<uint32_t, uint32_t>, BinData> bin_results;
  // Vector to store the formatted output rows before writing
  std::vector<BinOutputRow> output_rows;

  // build contig to intervals mapping for efficient lookup
  void build_contig_to_intervals_map();

  // process a single alignment into bin results (thread-safe)
  void process_single_alignment(const Alignment& aln, std::map<std::pair<uint32_t, uint32_t>, BinData>& target_bin_results);

  // merge thread-local BinData into global bin_results
  void merge_bin_data(const std::map<std::pair<uint32_t, uint32_t>, BinData>& local_data);

  // efficiently calculate position coverage only for positions with variants/clips
  void calc_position_coverage(std::map<std::pair<uint32_t, uint32_t>, BinData>& target_bin_results);

  // generate the output rows
  void generate_output_rows();

  public:
  QueryBin(
      const std::vector<Interval>& intervals,
      const AlignmentStore& store,
      int binsize,
      double seg_threshold = 0.2,
      double non_ref_threshold = 0.9,
      int num_threads = 0,
      ClipMode clip_mode = ClipMode::ALL,
      int clip_margin = 10,
      double min_mutations_percent = 0.0,
      double max_mutations_percent = 10.0,
      int min_alignment_length = 0,
      int max_alignment_length = 0,
      int min_seg_support = 2);

  // execute the query
  void execute();
  
  // execute data aggregation
  void aggregate_data();

  // write the output rows to a table
  void write_to_csv(const std::string& odir);

  // Getter for R interface
  const std::vector<BinOutputRow>& get_output_rows() const { return output_rows; }
};

#endif // QUERYBIN_H