#include "alignment_store.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <vector>

using std::cerr;
using std::endl;
using std::ifstream;
using std::ios;
using std::ofstream;
using std::string;
using std::unordered_map;
using std::vector;

// Helper function to write string to binary file
static void write_string(std::ofstream& file, const string& str)
{
  size_t len = str.size();
  file.write(reinterpret_cast<const char*>(&len), sizeof(len));
  file.write(str.c_str(), len);
}

// Helper function to read string from binary file
static string read_string(std::ifstream& file)
{
  size_t len;
  file.read(reinterpret_cast<char*>(&len), sizeof(len));
  string str(len, '\0');
  file.read(&str[0], len);
  return str;
}

const Mutation& AlignmentStore::get_mutation(uint32_t contig_idx, uint32_t mutation_idx) const
{
  auto contig_it = mutations_.find(contig_idx);
  massert(contig_it != mutations_.end(), "contig index %u not found in mutation store", contig_idx);
  massert(mutation_idx < contig_it->second.size(), "mutation index %u out of bounds for contig %u (size %zu)",
      mutation_idx, contig_idx, contig_it->second.size());
  return contig_it->second[mutation_idx];
}

void AlignmentStore::save(const string& filename)
{
  ofstream file(filename, ios::binary);
  massert(file.is_open(), "error opening file for writing: %s", filename.c_str());

  // Define a magic number/version for the format
  const string MAGIC_NUMBER = "ALNSTV2";
  file.write(MAGIC_NUMBER.c_str(), MAGIC_NUMBER.size());

  // Save contigs
  size_t num_contigs = contigs_.size();
  file.write(reinterpret_cast<const char*>(&num_contigs), sizeof(num_contigs));
  for (const auto& contig : contigs_) {
    write_string(file, contig.id);
    file.write(reinterpret_cast<const char*>(&contig.length), sizeof(contig.length));
  }

  // Save reads
  size_t num_reads = reads_.size();
  file.write(reinterpret_cast<const char*>(&num_reads), sizeof(num_reads));
  for (const auto& read : reads_) {
    write_string(file, read.id);
    file.write(reinterpret_cast<const char*>(&read.length), sizeof(read.length));
  }

  // Save mutations_ map
  size_t num_contigs_with_mutations = mutations_.size();
  file.write(reinterpret_cast<const char*>(&num_contigs_with_mutations), sizeof(num_contigs_with_mutations));

  for (const auto& pair : mutations_) {
    uint32_t contig_index = pair.first;
    const auto& mutations_vec = pair.second;
    size_t num_mutations_for_contig = mutations_vec.size();

    file.write(reinterpret_cast<const char*>(&contig_index), sizeof(contig_index));
    file.write(reinterpret_cast<const char*>(&num_mutations_for_contig), sizeof(num_mutations_for_contig));

    // Save each mutation for this contig
    for (const auto& mutation : mutations_vec) {
      // Write type (as underlying type, likely char or int8)
      file.write(reinterpret_cast<const char*>(&mutation.type), sizeof(mutation.type));
      // Always save position as uint32_t
      file.write(reinterpret_cast<const char*>(&mutation.position), sizeof(mutation.position));
      // Write the single nts string
      write_string(file, mutation.nts);
    }
  }

  // Save alignments
  size_t num_alignments = alignments_.size();
  file.write(reinterpret_cast<const char*>(&num_alignments), sizeof(num_alignments));
  for (const auto& alignment : alignments_) {
    // Write basic alignment data
    file.write(reinterpret_cast<const char*>(&alignment.read_index), sizeof(alignment.read_index));
    file.write(reinterpret_cast<const char*>(&alignment.contig_index), sizeof(alignment.contig_index));
    file.write(reinterpret_cast<const char*>(&alignment.read_start), sizeof(alignment.read_start));
    file.write(reinterpret_cast<const char*>(&alignment.read_end), sizeof(alignment.read_end));
    file.write(reinterpret_cast<const char*>(&alignment.contig_start), sizeof(alignment.contig_start));
    file.write(reinterpret_cast<const char*>(&alignment.contig_end), sizeof(alignment.contig_end));
    file.write(reinterpret_cast<const char*>(&alignment.is_reverse), sizeof(alignment.is_reverse));

    // Write mutation indices
    size_t num_mutation_indices = alignment.mutations.size(); // Now vector<uint32_t>
    file.write(reinterpret_cast<const char*>(&num_mutation_indices), sizeof(num_mutation_indices));
    for (uint32_t mutation_index : alignment.mutations) {
      file.write(reinterpret_cast<const char*>(&mutation_index), sizeof(mutation_index));
    }
  }

  file.close();

  // Set loaded flag to prevent further mutation additions via add_mutation
  loaded_ = true;

  // Organize alignments after loading
  organize_alignments();
}

void AlignmentStore::load(const string& filename)
{
  ifstream file(filename, ios::binary);
  massert(file.is_open(), "error opening file for reading: %s", filename.c_str());

  // Verify magic number
  const string EXPECTED_MAGIC = "ALNSTV2";
  char magic_buffer[8];
  file.read(magic_buffer, EXPECTED_MAGIC.size());
  massert(file.good() && string(magic_buffer, EXPECTED_MAGIC.size()) == EXPECTED_MAGIC,
      "invalid file format or version: %s", filename.c_str());

  // Clear existing data
  contigs_.clear();
  reads_.clear();
  alignments_.clear();
  mutations_.clear(); // Clear the new mutation store
  mutation_key_to_index_.clear(); // Clear the transient lookup map
  read_id_to_index.clear();
  contig_id_to_index.clear();
  alignment_index_by_contig_.clear();
  max_alignment_length_ = 0;

  // Load contigs
  size_t num_contigs;
  file.read(reinterpret_cast<char*>(&num_contigs), sizeof(num_contigs));
  contigs_.reserve(num_contigs);
  for (size_t i = 0; i < num_contigs; ++i) {
    string id = read_string(file);
    uint32_t length;
    file.read(reinterpret_cast<char*>(&length), sizeof(length));
    contigs_.emplace_back(id, length);
    contig_id_to_index[id] = i;
  }

  // Load reads
  size_t num_reads;
  file.read(reinterpret_cast<char*>(&num_reads), sizeof(num_reads));
  reads_.reserve(num_reads);
  for (size_t i = 0; i < num_reads; ++i) {
    string id = read_string(file);
    uint32_t length;
    file.read(reinterpret_cast<char*>(&length), sizeof(length));
    reads_.emplace_back(id, length);
    read_id_to_index[id] = i;
  }

  // Load mutations_ map
  size_t num_contigs_with_mutations;
  file.read(reinterpret_cast<char*>(&num_contigs_with_mutations), sizeof(num_contigs_with_mutations));

  for (size_t i = 0; i < num_contigs_with_mutations; ++i) {
    uint32_t contig_index;
    file.read(reinterpret_cast<char*>(&contig_index), sizeof(contig_index));

    size_t num_mutations_for_contig;
    file.read(reinterpret_cast<char*>(&num_mutations_for_contig), sizeof(num_mutations_for_contig));

    // Prepare vector for this contig's mutations
    vector<Mutation> mutations_vec;
    mutations_vec.reserve(num_mutations_for_contig);

    for (size_t j = 0; j < num_mutations_for_contig; ++j) {
      MutationType type;
      file.read(reinterpret_cast<char*>(&type), sizeof(type));

      // Always load position as uint32_t
      uint32_t position;
      file.read(reinterpret_cast<char*>(&position), sizeof(position));

      // Read nts string
      string nts = read_string(file);

      // Create mutation and add to vector using the new constructor
      mutations_vec.emplace_back(type, position, nts);
    }
    // Insert the loaded mutations into the map
    mutations_[contig_index] = std::move(mutations_vec);
  }

  // Load alignments
  size_t num_alignments;
  file.read(reinterpret_cast<char*>(&num_alignments), sizeof(num_alignments));
  alignments_.reserve(num_alignments);
  for (size_t i = 0; i < num_alignments; ++i) {
    Alignment alignment;

    // Read basic alignment data
    file.read(reinterpret_cast<char*>(&alignment.read_index), sizeof(alignment.read_index));
    file.read(reinterpret_cast<char*>(&alignment.contig_index), sizeof(alignment.contig_index));
    file.read(reinterpret_cast<char*>(&alignment.read_start), sizeof(alignment.read_start));
    file.read(reinterpret_cast<char*>(&alignment.read_end), sizeof(alignment.read_end));
    file.read(reinterpret_cast<char*>(&alignment.contig_start), sizeof(alignment.contig_start));
    file.read(reinterpret_cast<char*>(&alignment.contig_end), sizeof(alignment.contig_end));
    file.read(reinterpret_cast<char*>(&alignment.is_reverse), sizeof(alignment.is_reverse));

    // Read mutation indices
    size_t num_mutation_indices;
    file.read(reinterpret_cast<char*>(&num_mutation_indices), sizeof(num_mutation_indices));
    alignment.mutations.reserve(num_mutation_indices); // Reserve space in vector<uint32_t>
    for (size_t j = 0; j < num_mutation_indices; ++j) {
      uint32_t mutation_index;
      file.read(reinterpret_cast<char*>(&mutation_index), sizeof(mutation_index));
      alignment.mutations.push_back(mutation_index);
    }

    alignments_.push_back(std::move(alignment)); // Use move constructor
  }

  file.close();

  // Set loaded flag to prevent further mutation additions via add_mutation
  loaded_ = true;

  // Organize alignments after loading
  organize_alignments();
}

void AlignmentStore::organize_alignments()
{
  alignment_index_by_contig_.clear();
  max_alignment_length_ = 0; // Initialize max length calculation

  // Pre-allocate buckets for each contig
  for (size_t i = 0; i < contigs_.size(); ++i) {
    alignment_index_by_contig_[i] = {};
  }

  for (size_t i = 0; i < alignments_.size(); ++i) {
    const auto& alignment = alignments_[i];

    // Add index to the map
    auto it = alignment_index_by_contig_.find(alignment.contig_index);
    massert(it != alignment_index_by_contig_.end(), "alignment references unknown contig index %zu", alignment.contig_index);
    it->second.push_back(i);

    // Calculate and update max alignment length
    massert(alignment.contig_end >= alignment.contig_start, "alignment with end < start found (index %zu)", i);
    uint32_t current_length = alignment.contig_end - alignment.contig_start; // Length is end - start
    if (current_length > max_alignment_length_) {
      max_alignment_length_ = current_length;
    }
  }

  // Sort the alignment indices within each contig's vector based on start position
  for (auto& pair : alignment_index_by_contig_) {
    auto& indices = pair.second;
    std::sort(indices.begin(), indices.end(),
        [this](size_t index_a, size_t index_b) {
          return alignments_[index_a].contig_start < alignments_[index_b].contig_start;
        });
  }

  std::cout << "max alignment length found: " << max_alignment_length_ << std::endl;
}

void AlignmentStore::export_tab_delimited(const string& prefix)
{
  // Write main alignments file
  string alignments_file = prefix + "_alignments.txt";
  std::cout << "writing alignments to: " << alignments_file << std::endl;
  ofstream alignments_out(alignments_file);
  massert(alignments_out.is_open(), "Failed to open file for writing: %s", alignments_file.c_str());

  // Write mutations file
  string mutations_file = prefix + "_mutations.txt";
  std::cout << "writing mutations to: " << mutations_file << std::endl;
  ofstream mutations_out(mutations_file);
  massert(mutations_out.is_open(), "Failed to open file for writing: %s", mutations_file.c_str());

  // Write headers
  alignments_out << "read_id\tread_start\tread_end\tcontig_id\tcontig_start\tcontig_end\tmutation_count\tis_reverse\n";
  mutations_out << "read_id\tcontig_id\tmutation_type\tcontig_position\tnts\n";

  // Write alignments and mutations
  for (const auto& alignment : alignments_) {
    // Get read and contig IDs
    const string& read_id = get_read_id(alignment.read_index);
    const string& contig_id = get_contig_id(alignment.contig_index);

    // Write alignment data
    alignments_out << read_id << "\t"
                   << alignment.read_start << "\t"
                   << alignment.read_end << "\t"
                   << contig_id << "\t"
                   << alignment.contig_start << "\t"
                   << alignment.contig_end << "\t"
                   << alignment.mutations.size() << "\t"
                   << (alignment.is_reverse ? "true" : "false") << "\n";

    // Write detailed mutation data by fetching from store
    for (uint32_t mutation_index : alignment.mutations) { // Iterate indices
      // Get the actual mutation object
      const Mutation& mutation = get_mutation(alignment.contig_index, mutation_index);

      string mutation_type_str;
      switch (mutation.type) {
      case MutationType::SUBSTITUTION:
        mutation_type_str = "SUB";
        break;
      case MutationType::INSERTION:
        mutation_type_str = "INS";
        break;
      case MutationType::DELETION:
        mutation_type_str = "DEL";
        break;
      default:
        massert(false, "Unknown mutation type");
        break;
      }

      // Mutation position is now absolute contig position
      mutations_out << read_id << "\t"
                    << contig_id << "\t"
                    << mutation_type_str << "\t"
                    << mutation.position << "\t" // Absolute contig position
                    << mutation.nts << "\n";
    }
  }

  alignments_out.close();
  mutations_out.close();
}

size_t AlignmentStore::get_read_index(const string& read_id)
{
  auto it = read_id_to_index.find(read_id);
  massert(it != read_id_to_index.end(), "read not found: %s", read_id.c_str());
  return it->second;
}

size_t AlignmentStore::get_contig_index(const string& contig_id) const
{
  auto it = contig_id_to_index.find(contig_id);
  massert(it != contig_id_to_index.end(), "contig not found: %s", contig_id.c_str());
  return it->second;
}

size_t AlignmentStore::add_or_get_read_index(const string& read_id, uint32_t length)
{
  auto it = read_id_to_index.find(read_id);
  if (it != read_id_to_index.end()) {
    return it->second;
  } else {
    size_t new_index = reads_.size();
    reads_.emplace_back(read_id, length);
    read_id_to_index[read_id] = new_index;
    return new_index;
  }
}

size_t AlignmentStore::add_or_get_contig_index(const string& contig_id, uint32_t length)
{
  auto it = contig_id_to_index.find(contig_id);
  if (it != contig_id_to_index.end()) {
    return it->second;
  } else {
    size_t new_index = contigs_.size();
    contigs_.emplace_back(contig_id, length);
    contig_id_to_index[contig_id] = new_index;
    return new_index;
  }
}

bool AlignmentStore::has_read_id(const string& read_id) const
{
  return read_id_to_index.find(read_id) != read_id_to_index.end();
}

bool AlignmentStore::has_contig_id(const string& contig_id) const
{
  return contig_id_to_index.find(contig_id) != contig_id_to_index.end();
}

const string& AlignmentStore::get_read_id(size_t read_index) const
{
  massert(read_index < reads_.size(), "read index out of bounds: %zu", read_index);
  return reads_[read_index].id;
}

const string& AlignmentStore::get_contig_id(size_t contig_index) const
{
  massert(contig_index < contigs_.size(), "contig index out of bounds: %zu", contig_index);
  return contigs_[contig_index].id;
}

std::vector<std::reference_wrapper<const Alignment>> AlignmentStore::get_alignments_intersecting_interval(const Interval& interval, int max_alignments) const
{
  std::vector<std::reference_wrapper<const Alignment>> result;

  auto contig_map_it = contig_id_to_index.find(interval.contig);
  if (contig_map_it == contig_id_to_index.end()) {
    return result;
  }
  size_t contig_index = contig_map_it->second;

  auto align_map_it = alignment_index_by_contig_.find(contig_index);
  // It's possible a contig exists but has no alignments, so don't assert here.
  if (align_map_it == alignment_index_by_contig_.end()) {
    return result; // No alignments for this contig
  }

  const auto& indices = align_map_it->second;
  if (indices.empty()) {
    return result;
  }

  uint32_t min_possible_start = (interval.start >= max_alignment_length_) ? (interval.start - max_alignment_length_ + 1) : 0;

  auto it_start = std::lower_bound(indices.begin(), indices.end(), min_possible_start,
      [this](size_t index, uint32_t min_start) {
        return alignments_[index].contig_start < min_start;
      });

  auto it_end = std::upper_bound(indices.begin(), indices.end(), interval.end,
      [this](uint32_t query_end, size_t index) {
        return query_end < alignments_[index].contig_start;
      });

  // first pass: count valid alignments
  size_t valid_count = 0;
  for (auto it = it_start; it != it_end; ++it) {
    massert(*it < alignments_.size(), "alignment index out of bounds: %zu", *it);
    const auto& alignment = alignments_[*it];
    if (alignment.contig_end >= interval.start) {
      valid_count++;
    }
  }
  // print if we need to subsample
  if (max_alignments > 0 && valid_count > static_cast<size_t>(max_alignments)) {
    std::cout << "subsampling " << valid_count << " alignments to " << max_alignments << std::endl;
  }
  
  // second pass: collect alignments with subsampling if needed
  for (auto it = it_start; it != it_end; ++it) {
    massert(*it < alignments_.size(), "alignment index out of bounds: %zu", *it);
    const auto& alignment = alignments_[*it];
    if (alignment.contig_end >= interval.start) {
      if (max_alignments > 0 && valid_count > static_cast<size_t>(max_alignments)) {
        // subsample: skip based on random number
        if (rand() % valid_count < static_cast<size_t>(max_alignments)) {
          result.push_back(std::cref(alignment));
        }
      } else {
        result.push_back(std::cref(alignment));
      }
    }
  }

  return result;
}

// Add unique mutation (during build phase only)
uint32_t AlignmentStore::add_mutation(uint32_t contig_index, const Mutation& mutation)
{
  massert(!loaded_, "cannot add mutations after store has been loaded");

  string key = mutation.create_key(contig_index);
  auto it = mutation_key_to_index_.find(key);

  if (it != mutation_key_to_index_.end()) {
    // Found existing mutation, return its index
    return it->second;
  } else {
    // New mutation, add to store and map
    auto& contig_mutations = mutations_[contig_index]; // Get or create vector
    contig_mutations.push_back(mutation);
    uint32_t new_index = contig_mutations.size() - 1;
    mutation_key_to_index_[key] = new_index;
    return new_index;
  }
}



std::vector<AlignmentStore::BreakPosition> AlignmentStore::find_break_positions(
    uint32_t window_size, double p_threshold, uint32_t min_reads) const
{
  massert(window_size > 0, "window size must be positive");
  massert(p_threshold > 0.0 && p_threshold <= 1.0, "p threshold must be between 0 and 1");
  massert(min_reads > 0, "min reads must be positive");
  
  std::vector<BreakPosition> all_results;
  
  // process each contig
  for (size_t contig_idx = 0; contig_idx < contigs_.size(); ++contig_idx) {
    const Contig& contig = contigs_[contig_idx];
    
    // skip empty contigs
    if (contig.length == 0) continue;
    
    // build event count vectors
    std::vector<uint32_t> left_events(contig.length, 0);
    std::vector<uint32_t> right_events(contig.length, 0);
    
    // check if contig has alignments
    auto align_it = alignment_index_by_contig_.find(contig_idx);
    if (align_it != alignment_index_by_contig_.end()) {
      // count events from alignments
      for (size_t align_idx : align_it->second) {
        const Alignment& aln = alignments_[align_idx];
        const Read& read = reads_[aln.read_index];
        
        // count left break events (reads coming from the left)
        // forward reads starting at position OR reverse reads ending at position
        if ((aln.read_start == 0 && !aln.is_reverse) || 
            (aln.read_end == read.length && aln.is_reverse)) {
          uint32_t pos = aln.is_reverse ? aln.contig_start : aln.contig_start;
          if (pos < contig.length) {
            left_events[pos]++;
          }
        }
        
        // count right break events (reads coming from the right)
        // forward reads ending at position OR reverse reads starting at position
        if ((aln.read_end == read.length && !aln.is_reverse) || 
            (aln.read_start == 0 && aln.is_reverse)) {
          uint32_t pos = aln.is_reverse ? aln.contig_end - 1 : aln.contig_end - 1;
          if (pos < contig.length) {
            right_events[pos]++;
          }
        }
      }
    }
    
    // compute prefix sums for efficient window queries
    std::vector<uint32_t> left_prefix(contig.length + 1, 0);
    std::vector<uint32_t> right_prefix(contig.length + 1, 0);
    for (uint32_t i = 0; i < contig.length; ++i) {
      left_prefix[i + 1] = left_prefix[i] + left_events[i];
      right_prefix[i + 1] = right_prefix[i] + right_events[i];
    }
    
    // helper lambda to test events for a given orientation
    auto test_events = [&](const std::vector<uint32_t>& events, 
                          const std::vector<uint32_t>& prefix, 
                          const string& orientation) {
      for (uint32_t pos = 0; pos < contig.length; ++pos) {
        uint32_t t = events[pos];
        if (t < min_reads) continue;
        
        // compute window bounds
        uint32_t half_window = window_size / 2;
        uint32_t w_start = (pos >= half_window) ? pos - half_window : 0;
        uint32_t w_end = std::min(pos + half_window + 1, contig.length);
        
        // adjust window to maintain size when possible
        uint32_t actual_size = w_end - w_start;
        if (actual_size < window_size && contig.length >= window_size) {
          if (w_start == 0) {
            w_end = std::min(window_size, contig.length);
          } else if (w_end == contig.length) {
            w_start = contig.length - window_size;
          }
          actual_size = w_end - w_start;
        }
        
        uint32_t n_window = prefix[w_end] - prefix[w_start];
        if (n_window == 0) continue;
        
        double p0 = 1.0 / actual_size;
        double e = (double)n_window / actual_size;
        double enrichment = e > 0.0 ? (double)t / e : 0.0;
        double pval = binomial_right_tail(n_window, p0, t);
        
        all_results.push_back({
          contig.id, pos + 1, orientation, t, e, enrichment, pval, 0.0
        });
      }
    };
    
    // test left and right break events
    test_events(left_events, left_prefix, "left");
    test_events(right_events, right_prefix, "right");
  }
  
  // apply bh correction
  apply_bh_correction(all_results,
    [](BreakPosition& bp) -> double& { return bp.pval; },
    [](BreakPosition& bp) -> double& { return bp.qval; });
  
  // filter by q-value threshold
  std::vector<BreakPosition> filtered_results;
  for (const auto& result : all_results) {
    if (result.qval <= p_threshold) {
      filtered_results.push_back(result);
    }
  }
  
  // sort by contig id, then by position
  std::sort(filtered_results.begin(), filtered_results.end(),
    [](const BreakPosition& a, const BreakPosition& b) {
      if (a.contig_id != b.contig_id) {
        return a.contig_id < b.contig_id;
      }
      return a.position < b.position;
    });
  
  return filtered_results;
}

void AlignmentStore::count_short_indels(int min_indel_length)
{
  // skip if already computed for this min_indel_length
  if (last_min_indel_length_ == min_indel_length) {
    return;
  }

  // traverse all alignments and count short indels
  for (auto& alignment : alignments_) {
    alignment.short_indel_count = 0;
    
    for (uint32_t mutation_index : alignment.mutations) {
      const Mutation& mutation = get_mutation(alignment.contig_index, mutation_index);
      
      // check if mutation is a short indel
      if ((mutation.type == MutationType::INSERTION || mutation.type == MutationType::DELETION) &&
          static_cast<int>(mutation.nts.length()) < min_indel_length) {
        alignment.short_indel_count++;
      }
    }
  }

  // cache the min_indel_length to avoid recomputation
  last_min_indel_length_ = min_indel_length;
}

bool AlignmentStore::is_short_indel(uint32_t contig_idx, uint32_t mutation_idx) const
{
  if (last_min_indel_length_ == -1) {
    // not computed yet, assume all indels are long
    return false;
  }

  const Mutation& mutation = get_mutation(contig_idx, mutation_idx);
  return (mutation.type == MutationType::INSERTION || mutation.type == MutationType::DELETION) &&
         static_cast<int>(mutation.nts.length()) < last_min_indel_length_;
}

void AlignmentStore::init_read_alignment_index()
{
  if (read_alignment_index_built_) {
    return; // already built
  }

  cout << "building read-to-alignments index for " << alignments_.size() << " alignments" << endl;
  
  // clear any existing data
  read_to_alignment_indices_.clear();
  
  // build map of read_index to vector of alignment indices
  for (size_t i = 0; i < alignments_.size(); ++i) {
    const Alignment& aln = alignments_[i];
    read_to_alignment_indices_[aln.read_index].push_back(i);
  }
  
  // sort alignment indices by read_start for each read
  for (auto& pair : read_to_alignment_indices_) {
    vector<size_t>& alignment_indices = pair.second;
    std::sort(alignment_indices.begin(), alignment_indices.end(),
              [this](size_t idx1, size_t idx2) {
                return alignments_[idx1].read_start < alignments_[idx2].read_start;
              });
  }
  
  read_alignment_index_built_ = true;
}

bool AlignmentStore::is_alignment_local(const Alignment& alignment, int clip_margin) const
{
  if (!read_alignment_index_built_) {
    cerr << "error: read alignment index not built. call init_read_alignment_index() first" << endl;
    exit(EXIT_FAILURE);
  }
  
  uint32_t read_index = alignment.read_index;
  uint32_t focal_contig_index = alignment.contig_index;
  
  // find all alignments for this read
  auto it = read_to_alignment_indices_.find(read_index);
  if (it == read_to_alignment_indices_.end() || it->second.empty()) {
    return false; // no alignments found for this read
  }
  
  const vector<size_t>& alignment_indices = it->second;
  const Read& read = reads_[read_index];
  
  // look for focal-contig alignments at start and end of read
  bool has_focal_at_start = false;
  bool has_focal_at_end = false;
  
  // loop over all alignments for this read
  for (size_t aln_idx : alignment_indices) {
    const Alignment& aln = alignments_[aln_idx];
    
    // only consider alignments on the focal contig (ignore other contigs)
    if (aln.contig_index != focal_contig_index) {
      continue;
    }
    
    // check if this focal-contig alignment starts near beginning of read
    if (aln.read_start <= static_cast<uint32_t>(clip_margin)) {
      has_focal_at_start = true;
    }
    
    // check if this focal-contig alignment ends near end of read
    if (aln.read_end >= (read.length - static_cast<uint32_t>(clip_margin))) {
      has_focal_at_end = true;
    }
    
    // early exit if both conditions are met
    if (has_focal_at_start && has_focal_at_end) {
      return true;
    }
  }
  
  // return true only if we found focal-contig alignments at both start and end
  return has_focal_at_start && has_focal_at_end;
}

uint32_t AlignmentStore::get_alignment_overlap(const Alignment& input_alignment) const
{
  if (!read_alignment_index_built_) {
    cerr << "error: read alignment index not built. call init_read_alignment_index() first" << endl;
    exit(EXIT_FAILURE);
  }
  
  uint32_t read_index = input_alignment.read_index;
  
  // find all alignments for this read
  auto it = read_to_alignment_indices_.find(read_index);
  if (it == read_to_alignment_indices_.end() || it->second.empty()) {
    return 0; // no alignments found for this read
  }
  
  const vector<size_t>& alignment_indices = it->second;
  if (alignment_indices.size() <= 1) {
    return 0; // no other alignments to overlap with
  }
  
  // use binary search to find potentially overlapping alignments
  // calculate earliest possible start position for overlapping alignments
  uint32_t min_possible_start = (input_alignment.read_start >= max_alignment_length_) ? 
      (input_alignment.read_start - max_alignment_length_ + 1) : 0;
  
  // find first alignment that could overlap
  auto it_start = std::lower_bound(alignment_indices.begin(), alignment_indices.end(), min_possible_start,
      [this](size_t index, uint32_t min_start) {
        return alignments_[index].read_start < min_start;
      });
  
  // find last alignment that could overlap (starts before input_alignment ends)
  auto it_end = std::upper_bound(alignment_indices.begin(), alignment_indices.end(), input_alignment.read_end,
      [this](uint32_t query_end, size_t index) {
        return query_end < alignments_[index].read_start;
      });
  
  uint32_t max_overlap = 0;
  
  // check overlap with alignments in the potential range
  for (auto it = it_start; it != it_end; ++it) {
    size_t aln_idx = *it;
    const Alignment& other_alignment = alignments_[aln_idx];
    
    // skip self
    if (&other_alignment == &input_alignment) {
      continue;
    }
    
    // calculate overlap in read coordinates
    uint32_t overlap_start = std::max(input_alignment.read_start, other_alignment.read_start);
    uint32_t overlap_end = std::min(input_alignment.read_end, other_alignment.read_end);
    
    // check if there is actual overlap
    if (overlap_end > overlap_start) {
      uint32_t overlap_length = overlap_end - overlap_start;
      if (overlap_length > max_overlap) {
        max_overlap = overlap_length;
      }
    }
  }
  
  return max_overlap;
}

const Alignment* AlignmentStore::get_alignment_in_read(uint32_t read_index, uint32_t start_pos, uint32_t end_pos) const
{
  if (!read_alignment_index_built_) {
    cerr << "error: read alignment index not built. call init_read_alignment_index() first" << endl;
    exit(EXIT_FAILURE);
  }
  
  // check if interval is valid
  if (start_pos >= end_pos) {
    return nullptr; // no space for alignment
  }
  
  // find all alignments for this read
  auto it = read_to_alignment_indices_.find(read_index);
  if (it == read_to_alignment_indices_.end() || it->second.empty()) {
    return nullptr; // no alignments found for this read
  }
  
  const vector<size_t>& alignment_indices = it->second;
  if (alignment_indices.empty()) {
    return nullptr;
  }
  
  // use binary search to find potentially overlapping alignments
  // calculate earliest possible start position for alignments that could fit in interval
  uint32_t min_possible_start = (start_pos >= max_alignment_length_) ? 
      (start_pos - max_alignment_length_ + 1) : 0;
  
  // find first alignment that could potentially fit in the interval
  auto it_start = std::lower_bound(alignment_indices.begin(), alignment_indices.end(), min_possible_start,
      [this](size_t index, uint32_t min_start) {
        return alignments_[index].read_start < min_start;
      });
  
  // find last alignment that starts before or at end_pos
  auto it_end = std::upper_bound(alignment_indices.begin(), alignment_indices.end(), end_pos,
      [this](uint32_t query_end, size_t index) {
        return query_end < alignments_[index].read_start;
      });
  
  const Alignment* best_alignment = nullptr;
  uint32_t best_length = 0;
  uint32_t best_mutation_count = UINT32_MAX;
  
  // check alignments in the potential range
  for (auto it = it_start; it != it_end; ++it) {
    size_t aln_idx = *it;
    const Alignment& aln = alignments_[aln_idx];
    
    // check if alignment is completely within the interval
    if (aln.read_start >= start_pos && aln.read_end <= end_pos) {
      uint32_t length = aln.read_end - aln.read_start;
      uint32_t mutation_count = aln.get_mutation_count();
      
      // select longest alignment, with fewest mutations as tiebreaker
      if (length > best_length || (length == best_length && mutation_count < best_mutation_count)) {
        best_alignment = &aln;
        best_length = length;
        best_mutation_count = mutation_count;
      }
    }
  }
  
  return best_alignment;
}
