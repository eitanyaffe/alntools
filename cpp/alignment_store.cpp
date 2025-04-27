#include "alignment_store.h"
#include "utils.h"
#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <vector>

using std::cerr;
using std::endl;
using std::ifstream;
using std::ios;
using std::ofstream;
using std::string;
using std::unordered_map;
using std::vector;

void AlignmentStore::save_position(std::ofstream& file, uint32_t position, MutationIntType type)
{
  if (type == MutationIntType::UINT16) {
    uint16_t pos16 = static_cast<uint16_t>(position);
    file.write(reinterpret_cast<const char*>(&pos16), sizeof(pos16));
  } else {
    file.write(reinterpret_cast<const char*>(&position), sizeof(position));
  }
}

uint32_t AlignmentStore::load_position(std::ifstream& file, MutationIntType type)
{
  if (type == MutationIntType::UINT16) {
    uint16_t pos16;
    file.read(reinterpret_cast<char*>(&pos16), sizeof(pos16));
    return static_cast<uint32_t>(pos16);
  } else {
    uint32_t position;
    file.read(reinterpret_cast<char*>(&position), sizeof(position));
    return position;
  }
}

void AlignmentStore::save(const string& filename)
{
  ofstream file(filename, ios::binary);
  if (!file.is_open()) {
    cerr << "error opening file for writing: " << filename << endl;
    abort();
  }

  // Save reads
  size_t num_reads = reads_.size();
  file.write(reinterpret_cast<const char*>(&num_reads), sizeof(num_reads));
  for (const auto& read : reads_) {
    // Write the read id length and id
    size_t id_length = read.id.size();
    file.write(reinterpret_cast<const char*>(&id_length), sizeof(id_length));
    file.write(read.id.c_str(), id_length);

    // Write the read length
    file.write(reinterpret_cast<const char*>(&read.length), sizeof(read.length));
  }

  // Save contigs
  size_t num_contigs = contigs_.size();
  file.write(reinterpret_cast<const char*>(&num_contigs), sizeof(num_contigs));
  for (const auto& contig : contigs_) {
    // Write the contig id length and id
    size_t id_length = contig.id.size();
    file.write(reinterpret_cast<const char*>(&id_length), sizeof(id_length));
    file.write(contig.id.c_str(), id_length);

    // Write the contig length
    file.write(reinterpret_cast<const char*>(&contig.length), sizeof(contig.length));
  }

  // Determine maximum mutation position
  uint32_t max_position = 0;
  for (const auto& alignment : alignments_) {
    for (const auto& mutation : alignment.mutations) {
      if (mutation.position > max_position) {
        max_position = mutation.position;
      }
    }
  }

  // Determine position integer type
  MutationIntType mutation_int_type = (max_position <= std::numeric_limits<uint16_t>::max())
      ? MutationIntType::UINT16
      : MutationIntType::UINT32;
  string type_str = (mutation_int_type == MutationIntType::UINT16) ? "UINT16" : "UINT32";
  cout << "mutation position integer type: " << type_str << endl;

  // Save the mutation int type
  file.write(reinterpret_cast<const char*>(&mutation_int_type), sizeof(mutation_int_type));

  // Save alignments
  size_t num_alignments = alignments_.size();
  file.write(reinterpret_cast<const char*>(&num_alignments), sizeof(num_alignments));
  for (const auto& alignment : alignments_) {
    // Write the alignment data
    file.write(reinterpret_cast<const char*>(&alignment.read_index), sizeof(alignment.read_index));
    file.write(reinterpret_cast<const char*>(&alignment.contig_index), sizeof(alignment.contig_index));
    file.write(reinterpret_cast<const char*>(&alignment.read_start), sizeof(alignment.read_start));
    file.write(reinterpret_cast<const char*>(&alignment.read_end), sizeof(alignment.read_end));

    // Write contig start and end
    file.write(reinterpret_cast<const char*>(&alignment.contig_start), sizeof(alignment.contig_start));
    file.write(reinterpret_cast<const char*>(&alignment.contig_end), sizeof(alignment.contig_end));

    // Write is_reverse
    file.write(reinterpret_cast<const char*>(&alignment.is_reverse), sizeof(alignment.is_reverse));

    // Write mutations
    size_t num_mutations = alignment.mutations.size();
    file.write(reinterpret_cast<const char*>(&num_mutations), sizeof(num_mutations));
    for (const auto& mutation : alignment.mutations) {
      file.write(reinterpret_cast<const char*>(&mutation.type), sizeof(mutation.type));

      // Save position using the determined integer type
      save_position(file, mutation.position, mutation_int_type);

      // Write read_nts
      size_t read_nts_length = mutation.read_nts.size();
      file.write(reinterpret_cast<const char*>(&read_nts_length), sizeof(read_nts_length));
      file.write(mutation.read_nts.c_str(), read_nts_length);

      // Write ref_nts
      size_t ref_nts_length = mutation.ref_nts.size();
      file.write(reinterpret_cast<const char*>(&ref_nts_length), sizeof(ref_nts_length));
      file.write(mutation.ref_nts.c_str(), ref_nts_length);
    }
  }

  file.close();
}

void AlignmentStore::load(const string& filename)
{
  ifstream file(filename, ios::binary);
  if (!file.is_open()) {
    cerr << "error opening file: " << filename << endl;
    abort();
  }

  // Clear existing data
  reads_.clear();
  contigs_.clear();
  alignments_.clear();
  read_id_to_index.clear();
  contig_id_to_index.clear();
  alignment_index_by_contig_.clear();
  max_alignment_length_ = 0; // Reset max alignment length

  // reads
  size_t num_reads;
  file.read(reinterpret_cast<char*>(&num_reads), sizeof(num_reads));
  reads_.reserve(num_reads);
  for (size_t i = 0; i < num_reads; ++i) {
    Read read;

    // Read the read id length and id
    size_t id_length;
    file.read(reinterpret_cast<char*>(&id_length), sizeof(id_length));
    read.id.resize(id_length);
    file.read(&read.id[0], id_length);

    // Read the read length
    file.read(reinterpret_cast<char*>(&read.length), sizeof(read.length));

    reads_.push_back(read);
    read_id_to_index[read.id] = i;
  }

  // Load contigs
  size_t num_contigs;
  file.read(reinterpret_cast<char*>(&num_contigs), sizeof(num_contigs));
  contigs_.reserve(num_contigs);
  for (size_t i = 0; i < num_contigs; ++i) {
    Contig contig;

    // Read the contig id length and id
    size_t id_length;
    file.read(reinterpret_cast<char*>(&id_length), sizeof(id_length));
    contig.id.resize(id_length);
    file.read(&contig.id[0], id_length);

    // Read the contig length
    file.read(reinterpret_cast<char*>(&contig.length), sizeof(contig.length));

    contigs_.push_back(contig);
    contig_id_to_index[contig.id] = i;
  }

  // Read mutation integer type
  MutationIntType mutation_int_type;
  file.read(reinterpret_cast<char*>(&mutation_int_type), sizeof(mutation_int_type));
  string type_str = (mutation_int_type == MutationIntType::UINT16) ? "UINT16" : "UINT32";
  cout << "mutation position integer type: " << type_str << endl;

  // Load alignments
  size_t num_alignments;
  file.read(reinterpret_cast<char*>(&num_alignments), sizeof(num_alignments));
  alignments_.reserve(num_alignments);
  for (size_t i = 0; i < num_alignments; ++i) {
    Alignment alignment;

    // Read the alignment data
    file.read(reinterpret_cast<char*>(&alignment.read_index), sizeof(alignment.read_index));
    file.read(reinterpret_cast<char*>(&alignment.contig_index), sizeof(alignment.contig_index));
    file.read(reinterpret_cast<char*>(&alignment.read_start), sizeof(alignment.read_start));
    file.read(reinterpret_cast<char*>(&alignment.read_end), sizeof(alignment.read_end));
    // Read contig start and end
    file.read(reinterpret_cast<char*>(&alignment.contig_start), sizeof(alignment.contig_start));
    file.read(reinterpret_cast<char*>(&alignment.contig_end), sizeof(alignment.contig_end));
    // Read is_reverse
    file.read(reinterpret_cast<char*>(&alignment.is_reverse), sizeof(alignment.is_reverse));

    // Read mutations
    size_t num_mutations;
    file.read(reinterpret_cast<char*>(&num_mutations), sizeof(num_mutations));
    alignment.mutations.reserve(num_mutations);
    for (size_t j = 0; j < num_mutations; ++j) {
      MutationType type;
      string read_nts;
      string ref_nts;

      file.read(reinterpret_cast<char*>(&type), sizeof(type));

      // Load position using the read integer type
      uint32_t position = load_position(file, mutation_int_type);

      // Read read_nts
      size_t read_nts_length;
      file.read(reinterpret_cast<char*>(&read_nts_length), sizeof(read_nts_length));
      read_nts.resize(read_nts_length);
      file.read(&read_nts[0], read_nts_length);

      // Read ref_nts
      size_t ref_nts_length;
      file.read(reinterpret_cast<char*>(&ref_nts_length), sizeof(ref_nts_length));
      ref_nts.resize(ref_nts_length);
      file.read(&ref_nts[0], ref_nts_length);

      // Create mutation with the read data
      Mutation mutation(type, position, read_nts, ref_nts);
      alignment.mutations.push_back(mutation);
    }

    alignments_.push_back(alignment);
  }

  file.close();

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
    // Assuming end is inclusive, length is end - start + 1
    // Check for potential underflow if end < start although that shouldn't happen for valid alignments
    massert(alignment.contig_end >= alignment.contig_start, "alignment with end < start found (index %zu)", i);
    uint32_t current_length = alignment.contig_end - alignment.contig_start + 1;
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
  mutations_out << "read_id\tcontig_id\tmutation_type\tcontig_position\tread_position\tread_nts\tref_nts\n";

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
    // Write detailed mutation data
    for (const auto& mutation : alignment.mutations) {
      string mutation_type;
      switch (mutation.type) {
      case MutationType::SUBSTITUTION:
        mutation_type = "SUB";
        break;
      case MutationType::INSERTION:
        mutation_type = "INS";
        break;
      case MutationType::DELETION:
        mutation_type = "DEL";
        break;
      default:
        massert(false, "Unknown mutation type");
        break;
      }

      // Calculate read position from contig position
      uint32_t read_position;
      if (alignment.is_reverse) {
        read_position = alignment.read_end - (mutation.position - alignment.contig_start);
      } else {
        read_position = alignment.read_start + (mutation.position - alignment.contig_start);
      }

      mutations_out << read_id << "\t"
                    << contig_id << "\t"
                    << mutation_type << "\t"
                    << mutation.position << "\t"
                    << read_position << "\t"
                    << mutation.read_nts << "\t"
                    << mutation.ref_nts << "\n";
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

size_t AlignmentStore::get_contig_index(const string& contig_id)
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

std::vector<std::reference_wrapper<const Alignment>> AlignmentStore::get_alignments_in_interval(const Interval& interval) const
{
  std::vector<std::reference_wrapper<const Alignment>> result;

  auto contig_map_it = contig_id_to_index.find(interval.contig);
  massert(contig_map_it != contig_id_to_index.end(), "contig not found: %s", interval.contig.c_str());
  size_t contig_index = contig_map_it->second;

  auto align_map_it = alignment_index_by_contig_.find(contig_index);
  massert(align_map_it != alignment_index_by_contig_.end(), "contig index not found: %zu", contig_index);

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

  for (auto it = it_start; it != it_end; ++it) {
    massert(*it < alignments_.size(), "alignment index out of bounds: %zu", *it);
    const auto& alignment = alignments_[*it];
    if (alignment.contig_end >= interval.start) {
      result.push_back(std::cref(alignment));
    }
  }

  return result;
}
