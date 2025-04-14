#pragma once

#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include "aln_types.h"
#include "alignment_store.h"

class PafReader {
private:
    void split_line(const std::string& line, char delimiter, std::vector<std::string>& fields);
    void parse_cs_string(const std::string& cs_string, Alignment& alignment);
    void compare_cigar_and_cs(const std::string& cigar, const std::string& cs, const std::string& read_id);
    
    std::string filename_;
    std::fstream file_;

public:
    PafReader(const std::string& filename);
    std::unique_ptr<Alignment> next();
    bool read_paf(const std::string& filename, AlignmentStore& store);
}; 