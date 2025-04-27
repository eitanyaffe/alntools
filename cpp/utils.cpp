#include "utils.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "aln_types.h"

using namespace std;

void massert(bool cond, const char* fmt, ...)
{
  if (cond) {
    return;
  }

  fprintf(stderr, "Error: ");

  va_list argp;
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);

  fprintf(stderr, "\n");
  exit(-1);
}

void mexit(const char* fmt, ...)
{
  fprintf(stderr, "Error: ");

  va_list argp;
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);

  fprintf(stderr, "\n");
  exit(-1);
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
string apply_mutations(const string& seq, const vector<Mutation>& mutations,
    const string& read_id, const string& contig_id)
{
  string result;
  size_t prev_pos = 0, current_pos = 0;
  bool error_found = false;

  cout << "applying " << mutations.size() << " mutations to read " << read_id
       << " on contig " << contig_id << endl;

  int count = 0;
  // Process each mutation in order
  for (const auto& mutation : mutations) {
    current_pos = mutation.position;
    count++;

    // Verify position is within bounds
    massert(current_pos <= seq.size(),
        "mutation %u position %u is outside fragment bounds %u for read "
        "%s, contig %s",
        count, current_pos, seq.size(), read_id.c_str(), contig_id.c_str());

    // Copy unchanged sequence up to this mutation
    size_t gap = current_pos - prev_pos;
    if (gap > 0)
      result.append(seq.substr(prev_pos, gap));

    string ref_nts = to_upper(mutation.ref_nts);
    string read_nts = to_upper(mutation.read_nts);

    // Apply mutation based on type
    switch (mutation.type) {
    case MutationType::SUBSTITUTION:

      // Verify reference bases match expected
      if (seq.substr(current_pos, ref_nts.size()) != ref_nts) {
        printf(
            "reference bases at position %zu do not match expected. "
            "expected: %s, found: %s\n",
            current_pos, ref_nts.c_str(),
            seq.substr(current_pos, ref_nts.size()).c_str());
        error_found = true;
      }
      // Verify substitution lengths match
      if (read_nts.size() != ref_nts.size()) {
        printf(
            "substitution lengths do not match at position %zu. read: %zu, "
            "ref: %zu\n",
            current_pos, read_nts.size(), ref_nts.size());
        error_found = true;
      }
      // Add the substituted bases
      result.append(read_nts);
      current_pos += read_nts.size();
      break;

    case MutationType::INSERTION:
      // Add inserted bases
      result.append(read_nts);
      break;

    case MutationType::DELETION:
      // Verify reference bases match expected
      string obs_nts = seq.substr(current_pos, ref_nts.size());
      if (obs_nts != ref_nts) {
        printf(
            "reference bases at position %zu do not match expected for "
            "deletion. expected: %s, found: %s\n",
            current_pos, obs_nts.c_str(), ref_nts.c_str());
        error_found = true;
      }
      current_pos += ref_nts.size();
      break;
    }
    prev_pos = current_pos;
  }

  // Append any remaining reference sequence
  if (current_pos < seq.size()) {
    result.append(seq.substr(current_pos));
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

string generate_cs_tag(const Alignment& alignment)
{
  string result;
  uint32_t current_pos = 0;

  for (size_t i = 0; i < alignment.mutations.size(); ++i) {
    const Mutation& mut = alignment.mutations[i];

    // add match segment if there's a gap
    uint32_t gap = mut.position - current_pos;
    if (gap > 0) {
      result += ":" + std::to_string(gap);
      current_pos += gap;
    }

    // add the mutation
    switch (mut.type) {
    case MutationType::SUBSTITUTION:
      result += "*" + to_lower(mut.ref_nts) + to_lower(mut.read_nts);
      current_pos = mut.position + 1; // substitution advances one position
      break;
    case MutationType::INSERTION:
      result += "+" + to_lower(mut.read_nts);
      // insertion doesn't advance position
      break;
    case MutationType::DELETION:
      result += "-" + to_lower(mut.ref_nts);
      current_pos = mut.position + mut.ref_nts.length();
      break;
    default:
      massert(false, "unknown mutation type");
      break;
    }
  }

  // add a final match segment if needed (assuming we know the reference length)
  uint32_t gap = alignment.contig_end - alignment.contig_start - current_pos;
  if (gap > 0)
    result += ":" + std::to_string(gap);

  return result;
}