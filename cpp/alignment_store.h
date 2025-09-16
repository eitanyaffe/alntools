#pragma once

#include "aln_types.h"
#include <cstdint>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

using std::string;
using std::unordered_map;
using std::vector;

class AlignmentStore {
  private:
  std::vector<Contig> contigs_;
  std::vector<Read> reads_;
  std::vector<Alignment> alignments_;
  std::map<uint32_t, std::vector<Mutation>> mutations_;
  unordered_map<string, size_t> read_id_to_index;
  unordered_map<string, size_t> contig_id_to_index;

  // Transient map for mutation deduplication during initial build
  std::map<string, uint32_t> mutation_key_to_index_;
  unordered_map<size_t, vector<size_t>> alignment_index_by_contig_;
  // read index to vector of alignment indices (sorted by read_start)
  unordered_map<uint32_t, vector<size_t>> read_to_alignment_indices_;
  uint32_t max_alignment_length_ = 0;
  bool loaded_ = false; // Flag to prevent additions after loading
  int last_min_indel_length_ = -1; // cache for count_short_indels optimization
  bool read_alignment_index_built_ = false; // flag to track if read index is built

  public:
  // Add methods
  void add_contig(const Contig& contig) { contigs_.push_back(contig); }
  void add_read(const Read& read) { reads_.push_back(read); }
  // Adds a unique mutation (handling deduplication) and returns its index.
  // Only usable before load() is called.
  uint32_t add_mutation(uint32_t contig_index, const Mutation& mutation);
  void add_alignment(const Alignment& alignment) { alignments_.push_back(alignment); }

  // Getter methods
  const std::vector<Contig>& get_contigs() const { return contigs_; }
  const std::vector<Read>& get_reads() const { return reads_; }
  const std::vector<Alignment>& get_alignments() const { return alignments_; }

  // Non-const getters for modification
  std::vector<Contig>& get_contigs() { return contigs_; }
  std::vector<Read>& get_reads() { return reads_; }
  std::vector<Alignment>& get_alignments() { return alignments_; }

  // Get a specific mutation object by its contig index and mutation index
  const Mutation& get_mutation(uint32_t contig_idx, uint32_t mutation_idx) const;

  void export_tab_delimited(const string& prefix);

  // Save and load methods
  void save(const string& filename);
  void load(const string& filename);

  // Organize alignments
  void organize_alignments();

  // Getter methods
  size_t get_alignment_count() const { return alignments_.size(); }
  size_t get_read_count() const { return reads_.size(); }

  // Add or get read index
  size_t add_or_get_read_index(const string& read_id, uint32_t length);
  size_t add_or_get_contig_index(const string& contig_id, uint32_t length);

  bool has_read_id(const string& read_id) const;
  bool has_contig_id(const string& contig_id) const;

  // Get read index
  size_t get_read_index(const string& read_id);
  size_t get_contig_index(const string& contig_id) const;

  // Get id by index
  const string& get_read_id(size_t read_index) const;
  const string& get_contig_id(size_t contig_index) const;

  // get alignments in a specific interval
  std::vector<std::reference_wrapper<const Alignment>> get_alignments_in_interval(const Interval& interval, int max_alignments = 0) const;

  // break position detection result
  struct BreakPosition {
    string contig_id;
    uint32_t position; // 1-based
    string orientation; // "left" or "right"
    uint32_t t; // events at position
    double e; // expected events per position
    double enrichment; // observed/expected ratio
    double pval; // raw p-value
    double qval; // bh-adjusted q-value
  };

  // find positions with excessive read start/end events
  std::vector<BreakPosition> find_break_positions(uint32_t window_size, double p_threshold, uint32_t min_reads = 1) const;

  // count short indels in all alignments and update short_indel_count
  void count_short_indels(int min_indel_length);

  // check if a mutation is a short indel based on current cached min_indel_length
  bool is_short_indel(uint32_t contig_idx, uint32_t mutation_idx) const;

  // build read-to-alignments index map (sorted by read_start)
  void init_read_alignment_index();

  // check if alignment is locally aligned (first/last alignments on read are on same contig)
  bool is_alignment_local(const Alignment& alignment, int clip_margin = 10) const;

  // accessor methods for read-to-alignment index
  bool is_read_alignment_index_built() const { return read_alignment_index_built_; }
  const std::unordered_map<uint32_t, std::vector<size_t>>& get_read_to_alignment_indices() const { return read_to_alignment_indices_; }
};
