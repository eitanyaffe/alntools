#ifndef QUERYFULL_H
#define QUERYFULL_H

#include "alignment_store.h"
#include "aln_types.h"
#include <cstdint>
#include <string>
#include <vector>

// height calculation style
enum class HeightStyle {
  BY_COORD_LEFT, // minimal height without overlap, sort by start position
  BY_COORD_RIGHT, // minimal height without overlap, sort by end position
  BY_MUTATIONS // sort by mutation density
};

struct FullOutputAlignments {
  uint64_t alignment_index;
  std::string read_id;
  int read_length;
  std::string contig_id;
  int read_start;
  int read_end;
  int contig_start;
  int contig_end;
  bool is_reverse;
  std::string cs_tag;
  int num_mutations;
  int height;
};

struct FullOutputMutations {
  uint64_t alignment_index;
  std::string read_id;
  std::string contig_id;
  MutationType type;
  int position;
  std::string desc;
  int height;
};

struct FullOutputReads {
  std::string read_id;
  std::string contig_id;
  int read_length;
  int span_start;
  int span_end;
  int total_aligned_length;
  int num_alignments;
  int num_mutations;
  int height;
  bool read_reversed;
};

class QueryFull {
  private:
  std::vector<Interval> intervals;
  const AlignmentStore& store;
  HeightStyle height_style;
  int max_alignments;
  std::string alignment_filter;
  double min_mutations_density;

  std::vector<FullOutputAlignments> output_alignments;
  std::vector<FullOutputMutations> output_mutations;
  std::vector<FullOutputReads> output_reads;

  void generate_output_data();

  // calculate heights for alignments based on selected style
  void calculate_heights();

  // helper methods for different height calculation styles
  void calculate_heights_by_coord(bool sort_by_start);
  void calculate_heights_by_mutations();

  // helper methods for read processing
  void collect_reads_from_alignments();
  void calculate_read_heights();
  void assign_alignment_heights_from_reads();

  // helper methods for binary search in mutation-based height calculation
  bool has_overlap(const std::vector<std::pair<int, int>>& intervals, int start, int end);
  void add_sorted_interval(std::vector<std::pair<int, int>>& intervals, int start, int end);

  public:
  QueryFull(const std::vector<Interval>& intervals,
      const AlignmentStore& store,
      HeightStyle height_style = HeightStyle::BY_COORD_LEFT,
      int max_alignments = 0,
      const std::string& alignment_filter = "all",
      double min_mutations_density = 0.0);

  // execute the query
  void execute();

  // write the output rows to a table
  void write_to_csv(const std::string& ofn_prefix);

  // getters
  const std::vector<FullOutputAlignments>& get_output_alignments() const;
  const std::vector<FullOutputMutations>& get_output_mutations() const;
  const std::vector<FullOutputReads>& get_output_reads() const;

  // set height calculation style
  void set_height_style(HeightStyle style);
  HeightStyle get_height_style() const;
};

#endif // QUERYFULL_H