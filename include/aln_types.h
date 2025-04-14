#ifndef ALN_TYPES_H
#define ALN_TYPES_H

#include <string>
#include <vector>
#include <cstdint>

// Basic data structures for alignment data
struct Contig {
    std::string id;
    uint32_t length;

    Contig() = default;
    Contig(const std::string& id, uint32_t length) : id(id), length(length) {}
};

struct Read {
    std::string id;
    uint32_t length;
    bool is_reverse;
    std::vector<size_t> alignment_indices;

    Read() = default;
    Read(const std::string& id, uint32_t length) : id(id), length(length), is_reverse(false) {}
};

struct Alignment {
    std::string query_name;
    std::string target_name;
    uint32_t query_length;
    uint32_t target_length;
    uint32_t query_start;
    uint32_t query_end;
    uint32_t target_start;
    uint32_t target_end;
    uint32_t alignment_length;
    uint32_t alignment_score;
    bool is_reverse;

    Alignment() = default;
    Alignment(const std::string& query_name, const std::string& target_name,
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
};

// Enum for strand orientation
enum class Strand {
    FORWARD,
    REVERSE
};

#endif // ALN_TYPES_H 