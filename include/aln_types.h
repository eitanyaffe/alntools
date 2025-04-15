#pragma once

#include <string>
#include <vector>
#include <cstdint>

using std::string;
using std::vector;

// Enum for mutation types
enum class MutationType {
    SUBSTITUTION,  // Base substitution
    INSERTION,     // Insertion of bases
    DELETION       // Deletion of bases
};

// Structure to represent a mutation
struct Mutation {
    MutationType type;     // Type of mutation
    uint32_t position;     // Position in the reference sequence
    string query_bases;    // Bases in the query sequence (for substitutions and insertions)
    string target_bases;   // Bases in the target sequence (for substitutions and deletions)
    
    Mutation(MutationType type, uint32_t position, const string& query_bases = "",
	     const string& target_bases = "")
        : type(type), position(position), query_bases(query_bases), target_bases(target_bases) {}
};

// Basic data structures for alignment data
struct Contig {
    string id;
    uint32_t length;

    Contig(const string& id = "", uint32_t length = 0)
      : id(id), length(length) {}
};

struct Read {
    string id;
    uint32_t length;
    
    Read(const string& id = "", uint32_t length = 0)
      : id(id), length(length) {}
};

struct Alignment {
    uint32_t read_index;
    uint32_t contig_index;
    uint32_t read_start;
    uint32_t read_end;
    uint32_t contig_start;
    uint32_t contig_end;
    bool is_reverse;
    vector<Mutation> mutations;  // Vector to store mutations
    
    Alignment(uint32_t read_idx = 0, uint32_t contig_idx = 0,
             uint32_t c_start = 0, uint32_t c_end = 0,
             uint32_t r_start = 0, uint32_t r_end = 0,
             bool is_rev = false)
        : read_index(read_idx), contig_index(contig_idx),
          read_start(r_start), read_end(r_end),
          contig_start(c_start), contig_end(c_end),
          is_reverse(is_rev) {}
    
    // Add a mutation to the alignment
    void add_mutation(const Mutation& mutation) {
        mutations.push_back(mutation);
    }
    
    // Clear all mutations
    void clear_mutations() {
        mutations.clear();
    }
};
