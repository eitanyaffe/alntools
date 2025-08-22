#pragma once

#include <algorithm>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "aln_types.h"

using namespace std;

void massert(bool cond, const char* fmt, ...);

void mexit(const char* fmt, ...);
string reverse_complement(std::string seq);

// Convert a string to uppercase
string to_upper(string str);
string to_lower(const string& input);

// read FASTA
void read_fasta(const string& filename,
    const unordered_set<string>& contig_ids,
    unordered_map<string, string>& contigs);

// read FASTQ
void read_fastq(const string& filename,
    const unordered_set<string>& read_ids,
    unordered_map<string, string>& reads);

// write FASTA
void write_fasta(const string& filename,
    unordered_map<string, string>& contigs);

// write FASTQ
void write_fastq(const string& filename,
    unordered_map<string, string>& reads);

class AlignmentStore; // Forward declaration

// Function to apply mutations to a contig fragment
string apply_mutations(const string& contig_fragment,
    const vector<uint32_t>& mutation_indices,
    const AlignmentStore& store,
    const Alignment& alignment,
    const string& read_id,
    const string& contig_id);

enum class FileType {
  FASTA,
  FASTQ,
  UNKNOWN
};

FileType get_file_type(const std::string& filename);

double get_file_size_mb(const std::string& filename);

void read_intervals(const std::string& filename, std::vector<Interval>& intervals);

string generate_cs_tag(const Alignment& alignment, const AlignmentStore& store);

// alignment filtering
enum class ClipMode {
  ALL,                  // allow all alignments
  COMPLETE,            // alignment covers all read from start to end
  ALLOW_ONE_SIDE_CLIP, // allow clipped on one side (start at read start or read end)
  ONLY_ONE_SIDE_CLIPPED, // show only alignments clipped on one side
  ONLY_TWO_SIDE_CLIPPED, // show only alignments clipped on both sides
  ONLY_CLIPPED         // show alignments clipped on one or both sides
};

ClipMode string_to_clip_mode(const std::string& mode);

bool passes_alignment_filter(const Alignment& alignment, 
                           const AlignmentStore& store,
                           ClipMode clip_mode,
                           int clip_margin,
                           double min_mutations_percent,
                           double max_mutations_percent,
                           int min_alignment_length = 0,
                           int max_alignment_length = 0);

// statistical functions
double binomial_right_tail(uint32_t n, double p, uint32_t k);

template<typename T, typename GetPval, typename GetQval>
void apply_bh_correction(std::vector<T>& results, 
                        GetPval get_pval,
                        GetQval get_qval)
{
  if (results.empty()) return;
  
  // sort by p-value
  std::sort(results.begin(), results.end(), 
    [&get_pval](T& a, T& b) {
      return get_pval(a) < get_pval(b);
    });
  
  size_t m = results.size();
  for (size_t i = 0; i < m; ++i) {
    get_qval(results[i]) = get_pval(results[i]) * m / (i + 1);
  }
  
  // ensure monotonicity
  for (int i = m - 2; i >= 0; --i) {
    get_qval(results[i]) = std::min(get_qval(results[i]), get_qval(results[i + 1]));
  }
}