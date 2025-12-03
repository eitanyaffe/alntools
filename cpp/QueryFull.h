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
  BREAK_ON_OVERLAP, // start new chunk if next alignment overlaps (next.start < current.end - max_margin)
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
  bool strand_flipped;  // true if same contig as first alignment in chunk AND different strand
  int64_t vstart;  // view coordinate start (clipped to interval)
  int64_t vend;    // view coordinate end (clipped to interval)
  bool clip_start; // true if vstart corresponds to a clipped position
  bool clip_end;   // true if vend corresponds to a clipped position
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

struct FullOutputChunk {
  std::string chunk_id;
  std::string read_id;
  int read_length;
  int total_aligned_length;
  int num_alignments;
  int num_mutations;
  int height;
  bool read_reversed;
  int64_t vstart;  // view coordinates for height calculation (spans contigs)
  int64_t vend;
};

class QueryFull : public QueryBase {
  private:
  const AlignmentStore& store;
  HeightStyle height_style;
  int max_alignments;
  int max_margin;
  ChunkType chunk_type;

  std::vector<FullOutputAlignments> output_alignments;
  std::vector<FullOutputMutations> output_mutations;
  std::vector<FullOutputChunk> output_chunks;
  std::vector<int> alignment_to_chunk_index;

  // map from contig_id to intervals for efficient vcoord lookup
  std::map<std::string, std::vector<const Interval*>> contig_id_to_intervals_map;

  void collect_alignments_from_intervals();
  void build_contig_id_to_intervals_map();


  // helper methods for chunk processing
  void collect_chunks_from_alignments();
  std::vector<std::vector<FullOutputAlignments*>> group_alignments_into_chunks(const std::vector<FullOutputAlignments*>& alignments);
  void build_chunk_objects(const std::string& read_id, const std::vector<std::vector<FullOutputAlignments*>>& chunks);
  void calculate_chunk_heights();
  void assign_alignment_heights_from_chunks();
  void calculate_chunk_heights_by_coord(bool sort_by_start);
  void calculate_chunk_heights_by_mutations();

  // helper methods for binary search in mutation-based height calculation
  bool has_overlap(const std::vector<std::pair<int64_t, int64_t>>& intervals, int64_t start, int64_t end);
  void add_sorted_interval(std::vector<std::pair<int64_t, int64_t>>& intervals, int64_t start, int64_t end);
  
  // clip alignment to interval and compute vcoords
  // returns false if alignment is completely outside interval (warning logged, skip alignment)
  bool clip_alignment_to_interval(const Alignment& aln, const Interval& interval,
                                  int64_t& vstart, int64_t& vend,
                                  bool& clip_start, bool& clip_end) const;

  // verification: check that chunks and alignments are properly matched
  void verify_chunks_and_alignments() const;

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
      int max_margin = 10,
      ChunkType chunk_type = ChunkType::BREAK_ON_OVERLAP);

  // execute the query
  void execute();

  // write outputs to directory (create fixed filenames inside odir)
  void write_to_csv(const std::string& odir);

  // getters
  const std::vector<FullOutputAlignments>& get_output_alignments() const;
  const std::vector<FullOutputMutations>& get_output_mutations() const;
  const std::vector<FullOutputChunk>& get_output_chunks() const;

  // set height calculation style
  void set_height_style(HeightStyle style);
  HeightStyle get_height_style() const;
};

#endif // QUERYFULL_H