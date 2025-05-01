#ifndef QUERYFULL_H
#define QUERYFULL_H

#include "alignment_store.h"
#include "aln_types.h"
#include <cstdint>
#include <string>
#include <vector>

// height calculation style
enum class HeightStyle {
  BY_COORD, // minimal height without overlap
  BY_MUTATIONS // sort by mutation density
};

struct FullOutputAlignments {
  uint64_t alignment_index;
  std::string read_id;
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

class QueryFull {
  private:
  std::vector<Interval> intervals;
  const AlignmentStore& store;
  HeightStyle height_style;

  std::vector<FullOutputAlignments> output_alignments;
  std::vector<FullOutputMutations> output_mutations;

  void generate_output_data();

  // calculate heights for alignments based on selected style
  void calculate_heights();

  // helper methods for different height calculation styles
  void calculate_heights_by_coord();
  void calculate_heights_by_mutations();

  // helper methods for binary search in mutation-based height calculation
  bool has_overlap(const std::vector<std::pair<int, int>>& intervals, int start, int end);
  void add_sorted_interval(std::vector<std::pair<int, int>>& intervals, int start, int end);

  public:
  QueryFull(const std::vector<Interval>& intervals,
      const AlignmentStore& store,
      HeightStyle height_style = HeightStyle::BY_COORD);

  // execute the query
  void execute();

  // write the output rows to a table
  void write_to_csv(const std::string& ofn_prefix);

  // getters
  const std::vector<FullOutputAlignments>& get_output_alignments() const;
  const std::vector<FullOutputMutations>& get_output_mutations() const;

  // set height calculation style
  void set_height_style(HeightStyle style);
  HeightStyle get_height_style() const;
};

#endif // QUERYFULL_H