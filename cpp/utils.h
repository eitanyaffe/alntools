#pragma once

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <unordered_set>

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