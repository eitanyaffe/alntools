#ifndef QUERYFULL_H
#define QUERYFULL_H

#include "QueryBase.h"
#include <cstdint>
#include <string>
#include <vector>

// height calculation style
enum class HeightStyle {
  BY_COORD_LEFT, // minimal height without overlap, sort by start position
  BY_COORD_RIGHT, // minimal height without overlap, sort by end position
  BY_MUTATIONS // sort by mutation density
};

// chunk definition style
enum class ChunkType {
  READ, // define chunk as entire read
  ALIGNMENT, // define one chunk per alignment
  BREAK_ON_OVERLAP, // start new chunk if next alignment overlaps (next.start < current.end - max_gap)
  BREAK_ON_GAP // break if next alignment has a large gap (as implemented now)
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
  std::string chunk_id;
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

struct FullOutputChunk {
  std::string chunk_id;
  std::string read_id;
  std::string contig_id;
  int read_length;
  int span_start; // contig coordinates for height calculation
  int span_end;   // contig coordinates for height calculation
  int total_aligned_length;
  int num_alignments;
  int num_mutations;
  int height;
  bool read_reversed;
};

class QueryFull : public QueryBase {
  private:
  const AlignmentStore& store;
  HeightStyle height_style;
  int max_alignments;
  int max_gap;
  ChunkType chunk_type;

  std::vector<FullOutputAlignments> output_alignments;
  std::vector<FullOutputMutations> output_mutations;
  std::vector<FullOutputReads> output_reads;
  std::vector<FullOutputChunk> output_chunks;
  std::vector<int> alignment_to_chunk_index;

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

  // helper methods for chunk processing
  void collect_chunks_from_alignments();
  void calculate_chunk_heights();
  void assign_alignment_heights_from_chunks();
  void calculate_chunk_heights_by_coord(bool sort_by_start);
  void calculate_chunk_heights_by_mutations();

  // helper methods for binary search in mutation-based height calculation
  bool has_overlap(const std::vector<std::pair<int, int>>& intervals, int start, int end);
  void add_sorted_interval(std::vector<std::pair<int, int>>& intervals, int start, int end);

  public:
  QueryFull(const std::vector<Interval>& intervals,
      const AlignmentStore& store,
      HeightStyle height_style = HeightStyle::BY_COORD_LEFT,
      int max_alignments = 0,
      ClipMode clip_mode = ClipMode::ALL,
      int clip_margin = 10,
      double min_mutations_percent = 0.0,
      double max_mutations_percent = 10.0,
      int min_alignment_length = 0,
      int max_alignment_length = 0,
      int max_gap = 10,
      ChunkType chunk_type = ChunkType::BREAK_ON_OVERLAP);

  // execute the query
  void execute();

  // write the output rows to a table
  void write_to_csv(const std::string& ofn_prefix);

  // getters
  const std::vector<FullOutputAlignments>& get_output_alignments() const;
  const std::vector<FullOutputMutations>& get_output_mutations() const;
  const std::vector<FullOutputReads>& get_output_reads() const;
  const std::vector<FullOutputChunk>& get_output_chunks() const;

  // set height calculation style
  void set_height_style(HeightStyle style);
  HeightStyle get_height_style() const;
};

#endif // QUERYFULL_H