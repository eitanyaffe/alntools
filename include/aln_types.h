#ifndef ALN_TYPES_H
#define ALN_TYPES_H

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
    
    Mutation(MutationType type, uint32_t position, const string& query_bases = "", const string& target_bases = "")
        : type(type), position(position), query_bases(query_bases), target_bases(target_bases) {}
};

// Basic data structures for alignment data
struct Contig {
    string id;
    uint32_t length;

    Contig() = default;
    Contig(const string& id, uint32_t length) : id(id), length(length) {}
};

struct Read {
    string id;
    uint32_t length;
    bool is_reverse;
    vector<size_t> alignment_indices;

    Read() = default;
    Read(const string& id, uint32_t length) : id(id), length(length), is_reverse(false) {}
};

struct Alignment {
    string query_name;
    string target_name;
    uint32_t query_length;
    uint32_t target_length;
    uint32_t query_start;
    uint32_t query_end;
    uint32_t target_start;
    uint32_t target_end;
    uint32_t alignment_length;
    uint32_t alignment_score;
    bool is_reverse;
    vector<Mutation> mutations;  // List of mutations in this alignment

    Alignment() = default;
    Alignment(const string& query_name, const string& target_name,
             uint32_t query_length, uint32_t target_length,
             uint32_t query_start, uint32_t query_end,
             uint32_t target_start, uint32_t target_end,
             uint32_t alignment_length, uint32_t alignment_score,
             bool is_reverse)
        : query_name(query_name), target_name(target_name),
          query_length(query_length), target_length(target_length),
          query_start(query_start), query_end(query_end),
          target_start(target_start), target_end(target_end),
          alignment_length(alignment_length), alignment_score(alignment_score),
          is_reverse(is_reverse) {}
    
    // Add a mutation to the alignment
    void add_mutation(const Mutation& mutation) {
        mutations.push_back(mutation);
    }
    
    // Clear all mutations
    void clear_mutations() {
        mutations.clear();
    }
};

// Enum for strand orientation
enum class Strand {
    FORWARD,
    REVERSE
};

#endif // ALN_TYPES_H 