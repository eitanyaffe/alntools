#pragma once

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "aln_types.h"

using std::string;
using std::unordered_map;
using std::unordered_set;

void massert(bool cond, const char* fmt, ...);
void mexit(const char *fmt, ...);

std::string reverse_complement(std::string seq);

// Convert a string to uppercase
string to_upper(string str);

// Function to read FASTA file and create a map of contig IDs to sequences
void read_fasta(const string &filename,
		const unordered_set<string> &contig_ids,
		unordered_map<string, string>& contigs);


// Function to read FASTQ file and create a map of read IDs to sequences
void read_fastq(const string &filename,
		const unordered_set<string> &read_ids,
		unordered_map<string, string>& reads);


// transform string to uppercase
string to_upper(string str);

// Function to apply mutations to a contig fragment
string apply_mutations(const string &contig_fragment,
                       const vector<Mutation> &mutations,
                       uint32_t contig_start, 
					   bool quit_on_error);
