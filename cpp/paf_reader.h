#pragma once

#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include "aln_types.h"
#include "alignment_store.h"

using std::string;
using std::vector;
using std::unique_ptr;
using std::fstream;

class PafReader {
private:
  unordered_map<string, string> m_reads;
  unordered_map<string, string> m_contigs;
  
  void split_line(const string& line, char delimiter, vector<string>& fields);
  void parse_cs_string(const string& cs_string, vector<char>& actions, vector<string>& values);
  bool add_mutations(const string& cs_string, Alignment& alignment);

  void verify_cs_string(const string& cs_string, Alignment& alignment, size_t line_number);

  // returns true if PAF alignment tag is correct (by applying it to contig and comparing to read)
  bool verify_alignment(Alignment& alignment, const string& read_id, const string& contig_id);
  
public:
  void read_paf(const string& filename, AlignmentStore& store, int max_reads, bool should_verify, bool quit_on_error);

  // optionally load reads and contigs, used for verification
  void load_reads_contigs(const string& ifn_reads,
			  const string& ifn_contigs);
}; 
