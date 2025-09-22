#include "utils.h"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "alignment_store.h"
#include "aln_types.h"

using namespace std;

void massert(bool cond, const char* fmt, ...)
{
  if (cond) {
    return;
  }

  va_list args1;
  va_start(args1, fmt);
  va_list args2;
  va_copy(args2, args1);
  int size = std::vsnprintf(nullptr, 0, fmt, args1) + 1;
  va_end(args1);

  std::vector<char> buffer(size);
  std::vsnprintf(buffer.data(), size, fmt, args2);
  va_end(args2);

  throw std::runtime_error("Assertion failed: " + std::string(buffer.data()));
}

void mexit(const char* fmt, ...)
{
  va_list args1;
  va_start(args1, fmt);
  va_list args2;
  va_copy(args2, args1);
  int size = std::vsnprintf(nullptr, 0, fmt, args1) + 1;
  va_end(args1);

  std::vector<char> buffer(size);
  std::vsnprintf(buffer.data(), size, fmt, args2);
  va_end(args2);

  throw std::runtime_error("Error: " + std::string(buffer.data()));
}

// Function to convert a string to uppercase
string to_upper(string str)
{
  std::transform(str.begin(), str.end(), str.begin(), ::toupper);
  return str;
}

string to_lower(const string& input)
{
  string result = input;
  transform(result.begin(), result.end(), result.begin(), ::tolower);
  return result;
}

std::string reverse_complement(std::string seq)
{
  std::string result = seq;
  int N = seq.length();
  for (int i = 0; i < N; i++) {
    char c = seq[N - i - 1], r;
    switch (c) {
    case 'A':
      r = 'T';
      break;
    case 'G':
      r = 'C';
      break;
    case 'C':
      r = 'G';
      break;
    case 'T':
      r = 'A';
      break;
    default:
      r = c;
    }
    result[i] = r;
  }
  return (result);
}

// Function to read FASTA file and create a map of contig IDs to sequences
void read_fasta(const string& filename, const unordered_set<string>& contig_ids,
    unordered_map<string, string>& contigs)
{
  cout << "Reading FASTA file: " << filename << endl;
  ifstream file(filename);
  string line, id, sequence;
  while (getline(file, line)) {
    if (line[0] == '>') {
      if (!id.empty() && (contig_ids.empty() || contig_ids.find(id) != contig_ids.end())) {
        contigs[id] = sequence;
      }
      id = line.substr(
          1, line.find_first_of(" ") - 1); // Extract ID before first space/tab
      sequence.clear();
    } else {
      sequence += line;
    }
  }
  if (!id.empty() && (contig_ids.empty() || contig_ids.find(id) != contig_ids.end())) {
    contigs[id] = sequence;
  }
  cout << "read " << contigs.size() << " contigs" << endl;
}

// Function to read FASTQ file and create a map of read IDs to sequences
void read_fastq(const string& filename, const unordered_set<string>& read_ids,
    unordered_map<string, string>& reads)
{
  cout << "Reading FASTQ file: " << filename << endl;
  ifstream file(filename);
  string line, id, sequence;
  while (getline(file, line)) {
    if (line[0] == '@') {
      id = line.substr(1);
      getline(file, sequence); // Read the sequence line
      if (read_ids.empty() || read_ids.find(id) != read_ids.end())
        reads[id] = sequence;
      getline(file, line); // Skip the '+' line
      getline(file, line); // Skip the quality line
    }
    if (!read_ids.empty() && read_ids.size() == reads.size())
      break;
  }
}

void write_fasta(const string& filename,
    unordered_map<string, string>& contigs)
{
  cout << "writing " << contigs.size() << " contigs to " << filename << endl;
  ofstream file(filename);
  if (!file.is_open()) {
    fprintf(stderr, "could not open file %s for writing\n", filename.c_str());
    exit(EXIT_FAILURE);
  }

  for (const auto& entry : contigs) {
    file << ">" << entry.first << "\n";
    file << entry.second << "\n";
  }

  file.close();
}

void write_fastq(const string& filename, unordered_map<string, string>& reads)
{
  cout << "writing " << reads.size() << " reads to " << filename << endl;
  ofstream file(filename);
  if (!file.is_open()) {
    fprintf(stderr, "could not open file %s for writing\n", filename.c_str());
    exit(EXIT_FAILURE);
  }

  for (const auto& entry : reads) {
    file << "@" << entry.first << "\n";
    file << entry.second << "\n";
    file << "+" << "\n";
    file << string(entry.second.length(), 'I')
         << "\n"; // Placeholder quality scores
  }

  file.close();
}

// Function to apply mutations to a contig fragment
// Note: Now takes AlignmentStore to fetch mutations by index
string apply_mutations(const string& seq, const vector<uint32_t>& mutation_indices,
    const AlignmentStore& store, const Alignment& alignment,
    const string& read_id, const string& contig_id)
{
  string result;
  size_t prev_pos_rel = 0, current_pos_rel = 0;
  bool error_found = false;

  // cout << "applying " << mutation_indices.size() << " mutations to read " << read_id
  //      << " on contig " << contig_id << " starting at " << alignment.contig_start << endl;

  int count = 0;
  // Process each mutation index
  for (uint32_t mut_idx : mutation_indices) {
    count++;
    // Fetch the actual mutation object
    const Mutation& mutation = store.get_mutation(alignment.contig_index, mut_idx);

    // Get the absolute position from the mutation object
    uint32_t current_pos_abs = mutation.position;

    // Calculate the position relative to the fragment start
    massert(current_pos_abs >= alignment.contig_start,
        "Mutation absolute position %u is before alignment start %u",
        current_pos_abs, alignment.contig_start);
    current_pos_rel = current_pos_abs - alignment.contig_start;

    // Verify relative position is within fragment bounds
    massert(current_pos_rel <= seq.size(),
        "mutation %d relative position %zu is outside fragment bounds %zu for read "
        "%s, contig %s (abs pos %u, aln start %u)",
        count, current_pos_rel, seq.size(), read_id.c_str(), contig_id.c_str(),
        current_pos_abs, alignment.contig_start);

    // Copy unchanged sequence up to this mutation
    size_t gap = current_pos_rel - prev_pos_rel;
    if (gap > 0) {
      massert(prev_pos_rel + gap <= seq.size(), "Gap calculation error: prev=%zu, gap=%zu, size=%zu", prev_pos_rel, gap, seq.size());
      result.append(seq.substr(prev_pos_rel, gap));
    }

    // Apply mutation based on type
    switch (mutation.type) {
    case MutationType::SUBSTITUTION: {
      massert(mutation.nts.length() == 2, "SUB mutation nts length is not 2");
      string read_nts = to_upper(mutation.nts.substr(0, 1)); // read is first char
      string ref_nts = to_upper(mutation.nts.substr(1, 1)); // ref is second char

      // Verify reference bases match expected
      massert(current_pos_rel + ref_nts.size() <= seq.size(), "Substitution check out of bounds");
      if (seq.substr(current_pos_rel, ref_nts.size()) != ref_nts) {
        printf(
            "reference bases at relative position %zu (abs %u) do not match expected. "
            "expected: %s, found: %s\n",
            current_pos_rel, current_pos_abs, ref_nts.c_str(),
            seq.substr(current_pos_rel, ref_nts.size()).c_str());
        error_found = true;
      }
      // Add the substituted bases (read_nts)
      result.append(read_nts);
      // Advance relative position past the substitution (length is always 1 for SUB)
      current_pos_rel += 1; // ref_nts.size() is always 1
      break;
    }
    case MutationType::INSERTION: {
      string read_nts = to_upper(mutation.nts);
      // Add inserted bases
      result.append(read_nts);
      // Position does not advance on reference
      break;
    }
    case MutationType::DELETION: {
      string ref_nts = to_upper(mutation.nts);
      // Verify reference bases match expected
      massert(current_pos_rel + ref_nts.size() <= seq.size(), "Deletion check out of bounds");
      string obs_nts = seq.substr(current_pos_rel, ref_nts.size());
      if (obs_nts != ref_nts) {
        printf(
            "reference bases at relative position %zu (abs %u) do not match expected for "
            "deletion. expected: %s, found: %s\n",
            current_pos_rel, current_pos_abs, obs_nts.c_str(), ref_nts.c_str());
        error_found = true;
      }
      // Advance relative position past the deletion
      current_pos_rel += ref_nts.size();
      break;
    }
    }
    // Update previous relative position for next iteration's gap calculation
    prev_pos_rel = current_pos_rel;
  }

  // Append any remaining reference sequence
  if (prev_pos_rel < seq.size()) {
    result.append(seq.substr(prev_pos_rel));
  }

  if (error_found) {
    mexit("error found, quitting", -1);
  }
  return result;
}

FileType get_file_type(const std::string& filename)
{
  cout << "getting file type for " << filename << endl;
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "failed to open file: " << filename << std::endl;
    return FileType::UNKNOWN;
  }

  std::string first_line;
  std::getline(file, first_line);
  file.close();

  if (first_line.empty()) {
    std::cerr << "file is empty: " << filename << std::endl;
    return FileType::UNKNOWN;
  }

  if (first_line[0] == '>') {
    return FileType::FASTA;
  } else if (first_line[0] == '@') {
    return FileType::FASTQ;
  } else {
    std::cerr << "unknown file format: " << filename << std::endl;
    return FileType::UNKNOWN;
  }
}

double get_file_size_mb(const std::string& filename)
{
  std::ifstream file(filename, std::ios::binary | std::ios::ate);
  if (!file.is_open()) {
    std::cerr << "failed to open file: " << filename << std::endl;
    return 0.0;
  }
  std::streamsize size = file.tellg();
  file.close();
  return static_cast<double>(size) / (1024.0 * 1024.0);
}

void read_intervals(const std::string& filename,
    std::vector<Interval>& intervals)
{
  std::ifstream file(filename);
  if (!file.is_open()) {
    cerr << "error: could not open file " << filename << " for reading" << endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  if (getline(file, line)) {
    // Verify header
    if (line != "contig\tstart\tend") {
      cerr << "error: invalid header in intervals file. Expected "
              "'contig\\tstart\\tend'"
           << endl;
      exit(EXIT_FAILURE);
    }
  }
  // read
  while (getline(file, line)) {
    std::istringstream iss(line);
    std::string contig;
    uint32_t start, end;
    if (!(iss >> contig >> start >> end)) {
      cerr << "error: malformed line in intervals file: " << line << endl;
      exit(EXIT_FAILURE);
    }
    intervals.emplace_back(contig, start, end);
  }
  file.close();
}

// Note: Needs update to work with mutation indices and AlignmentStore
string generate_cs_tag(const Alignment& alignment, const AlignmentStore& store)
{
  string result;
  // Position relative to the start of the alignment
  uint32_t current_relative_pos = 0;

  for (const auto& mut_idx : alignment.mutations) {
    // Fetch the actual mutation object
    const Mutation& mut = store.get_mutation(alignment.contig_index, mut_idx);

    // Calculate the relative position for this mutation
    massert(mut.position >= alignment.contig_start, "Mutation position %u before alignment start %u", mut.position, alignment.contig_start);
    uint32_t mutation_relative_pos = mut.position - alignment.contig_start;

    // add match segment if there's a gap relative to the previous operation
    uint32_t gap = mutation_relative_pos - current_relative_pos;
    if (gap > 0) {
      result += ":" + std::to_string(gap);
      current_relative_pos += gap;
    }

    // add the mutation
    switch (mut.type) {
    case MutationType::SUBSTITUTION: {
      massert(mut.nts.length() == 2, "SUB mutation nts length is not 2 for cs tag generation");
      string read_nt = mut.nts.substr(0, 1);
      string ref_nt = mut.nts.substr(1, 1);
      result += "*" + to_lower(ref_nt) + to_lower(read_nt);
      // Substitution advances relative position by 1
      current_relative_pos = mutation_relative_pos + 1;
      break;
    }
    case MutationType::INSERTION:
      result += "+" + to_lower(mut.nts);
      // insertion doesn't advance relative position
      break;
    case MutationType::DELETION: {
      string deleted_nts = mut.nts;
      result += "-" + to_lower(deleted_nts);
      // Deletion advances relative position by the length of the deleted sequence
      current_relative_pos = mutation_relative_pos + deleted_nts.length();
      break;
    }
    default:
      massert(false, "unknown mutation type");
      break;
    }
  }

  uint32_t gap = alignment.contig_end - alignment.contig_start - current_relative_pos;
  if (gap > 0)
    result += ":" + std::to_string(gap);

  return result;
}

// calculate binomial right-tail probability
double binomial_right_tail(uint32_t n, double p, uint32_t k)
{
  if (k > n) return 0.0;
  if (p <= 0.0) return (k == 0) ? 1.0 : 0.0;
  if (p >= 1.0) return (k <= n) ? 1.0 : 0.0;
  
  // use normal approximation for large n
  if (n > 100) {
    double mu = n * p;
    double sigma = sqrt(n * p * (1 - p));
    if (sigma > 0) {
      double z = (k - 0.5 - mu) / sigma; // continuity correction
      return 0.5 * erfc(z / sqrt(2.0));
    }
  }
  
  // exact calculation using log probabilities
  double log_sum = -std::numeric_limits<double>::infinity();
  for (uint32_t i = k; i <= n; ++i) {
    double log_prob = lgamma(n + 1) - lgamma(i + 1) - lgamma(n - i + 1) 
                     + i * log(p) + (n - i) * log(1 - p);
    if (log_sum == -std::numeric_limits<double>::infinity()) {
      log_sum = log_prob;
    } else {
      double max_log = std::max(log_sum, log_prob);
      log_sum = max_log + log1p(exp(std::min(log_sum, log_prob) - max_log));
    }
  }
  return exp(log_sum);
}

// alignment filtering functions
ClipMode string_to_clip_mode(const std::string& mode)
{
  if (mode == "all") {
    return ClipMode::ALL;
  } else if (mode == "complete") {
    return ClipMode::COMPLETE;
  } else if (mode == "allow_one_side_clip") {
    return ClipMode::ALLOW_ONE_SIDE_CLIP;
  } else if (mode == "only_one_side_clipped") {
    return ClipMode::ONLY_ONE_SIDE_CLIPPED;
  } else if (mode == "only_two_side_clipped") {
    return ClipMode::ONLY_TWO_SIDE_CLIPPED;
  } else if (mode == "only_clipped") {
    return ClipMode::ONLY_CLIPPED;
  } else if (mode == "local_align") {
    return ClipMode::LOCAL_ALIGN;
  } else {
    cerr << "warning: invalid clip_mode '" << mode << "', defaulting to 'all'" << endl;
    return ClipMode::ALL;
  }
}

bool passes_alignment_filter(const Alignment& alignment, 
                           const AlignmentStore& store,
                           ClipMode clip_mode,
                           int clip_margin,
                           double min_mutations_percent,
                           double max_mutations_percent,
                           int min_alignment_length,
                           int max_alignment_length)
{
  // get read length
  uint32_t read_length = store.get_reads()[alignment.read_index].length;
  
  // check clipping based on clip_mode
  bool starts_at_beginning = (alignment.read_start <= static_cast<uint32_t>(clip_margin));
  bool ends_at_end = (alignment.read_end >= (read_length - static_cast<uint32_t>(clip_margin)));
  
  if (clip_mode == ClipMode::COMPLETE) {
    // alignment must cover all read from start to end (with margin)
    if (!starts_at_beginning || !ends_at_end) {
      return false;
    }
  } else if (clip_mode == ClipMode::ALLOW_ONE_SIDE_CLIP) {
    // allow clipped on one side: either start at read start or end at read end
    if (!starts_at_beginning && !ends_at_end) {
      return false;
    }
  } else if (clip_mode == ClipMode::ONLY_ONE_SIDE_CLIPPED) {
    // show only alignments clipped on exactly one side
    if ((starts_at_beginning && ends_at_end) || (!starts_at_beginning && !ends_at_end)) {
      return false;
    }
  } else if (clip_mode == ClipMode::ONLY_TWO_SIDE_CLIPPED) {
    // show only alignments clipped on both sides
    if (starts_at_beginning || ends_at_end) {
      return false;
    }
  } else if (clip_mode == ClipMode::ONLY_CLIPPED) {
    // show alignments clipped on one or both sides
    if (starts_at_beginning && ends_at_end) {
      return false;
    }
  } else if (clip_mode == ClipMode::LOCAL_ALIGN) {
    // show only locally aligned reads (first/last alignments on same contig)
    if (!store.is_alignment_local(alignment, clip_margin)) {
      return false;
    }
  }
  // ClipMode::ALL allows all alignments, no clipping check needed
  
  // check mutations percentage range
  if (min_mutations_percent >= 0.0 || max_mutations_percent >= 0.0) {
    uint32_t alignment_length = alignment.contig_end - alignment.contig_start;
    if (alignment_length > 0) {
      double actual_mutations_percent = (static_cast<double>(alignment.get_mutation_count()) / alignment_length) * 100.0;
      
      // check minimum threshold
      if (min_mutations_percent >= 0.0 && actual_mutations_percent < min_mutations_percent) {
        return false;
      }
      
      // check maximum threshold
      if (max_mutations_percent >= 0.0) {
        // special case: if max_mutations_percent is 0, only allow alignments with exactly 0 mutations
        if (max_mutations_percent == 0.0) {
          if (alignment.get_mutation_count() > 0) {
            return false;
          }
        } else {
          // normal case: filter if above threshold
          if (actual_mutations_percent > max_mutations_percent) {
            return false;
          }
        }
      }
    }
  }
  
  // check alignment length (based on read coordinates)
  if (min_alignment_length > 0 || max_alignment_length > 0) {
    uint32_t alignment_length = alignment.read_end - alignment.read_start;
    
    // check minimum length
    if (min_alignment_length > 0 && alignment_length < static_cast<uint32_t>(min_alignment_length)) {
      return false;
    }
    
    // check maximum length (0 means no limit)
    if (max_alignment_length > 0 && alignment_length > static_cast<uint32_t>(max_alignment_length)) {
      return false;
    }
  }
  
  return true;
}

void init_local_align_if_needed(AlignmentStore& store, ClipMode clip_mode)
{
  if (clip_mode == ClipMode::LOCAL_ALIGN) {
    store.init_read_alignment_index();
  }
}

bool should_skip_short_indel(const AlignmentStore& store, uint32_t contig_index, uint32_t mutation_index)
{
  return store.is_short_indel(contig_index, mutation_index);
}

void safe_open_file_for_writing(const std::string& filename, std::ofstream& file)
{
  file.open(filename);
  if (!file.is_open()) {
    cerr << "error: cannot open file " << filename << " for writing" << endl;
    exit(1);
  }
}

std::map<std::string, std::string> read_library_table(const std::string& filename)
{
    std::map<std::string, std::string> library_files;
    
    std::ifstream lib_file(filename);
    if (!lib_file.is_open()) {
        std::cerr << "error: cannot open libraries file " << filename << std::endl;
        std::exit(1);
    }
    
    std::string line;
    bool first_line = true;
    while (std::getline(lib_file, line)) {
        if (first_line) {
            first_line = false;
            continue; // skip header
        }
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::string lib_id, aln_file;
        if (std::getline(iss, lib_id, '\t') && std::getline(iss, aln_file, '\t')) {
            library_files[lib_id] = aln_file;
        }
    }
    lib_file.close();
    
    std::cout << "found " << library_files.size() << " libraries in " << filename << std::endl;
    return library_files;
}

std::map<std::string, LibraryInfo> read_library_table_extended(const std::string& filename)
{
    std::map<std::string, LibraryInfo> library_info;
    
    std::ifstream lib_file(filename);
    if (!lib_file.is_open()) {
        std::cerr << "error: cannot open libraries file " << filename << std::endl;
        std::exit(1);
    }
    
    std::string line;
    bool first_line = true;
    while (std::getline(lib_file, line)) {
        if (first_line) {
            first_line = false;
            continue; // skip header
        }
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::string lib_id, aln_file, read_file;
        
        // read required columns
        if (std::getline(iss, lib_id, '\t') && std::getline(iss, aln_file, '\t')) {
            // read optional third column (read file)
            if (!std::getline(iss, read_file, '\t')) {
                read_file = "";  // no read file specified
            }
            
            library_info[lib_id] = LibraryInfo(lib_id, aln_file, read_file);
        }
    }
    lib_file.close();
    
    std::cout << "found " << library_info.size() << " libraries in " << filename << std::endl;
    
    // report which libraries have read files
    int with_reads = 0;
    for (const auto& entry : library_info) {
        if (!entry.second.read_file.empty()) {
            with_reads++;
        }
    }
    std::cout << "  " << with_reads << " libraries have read files specified" << std::endl;
    
    return library_info;
}

