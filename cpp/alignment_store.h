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
  uint32_t max_alignment_length_ = 0;
  bool loaded_ = false; // Flag to prevent additions after loading

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

  // Get read index
  size_t get_read_index(const string& read_id);
  size_t get_contig_index(const string& contig_id) const;

  // Get id by index
  const string& get_read_id(size_t read_index) const;
  const string& get_contig_id(size_t contig_index) const;

  // New method to get alignments in a specific interval
  std::vector<std::reference_wrapper<const Alignment>> get_alignments_in_interval(const Interval& interval) const;
};
