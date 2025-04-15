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
    void split_line(const string& line, char delimiter, vector<string>& fields);
    void parse_cs_string(const string& cs_string, vector<char>& actions, vector<string>& values);
    size_t add_mutations(const string& cs_string, Alignment& alignment);
  
public:
  void read_paf(const string& filename, AlignmentStore& store, int max_reads);
}; 
