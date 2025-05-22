#pragma once

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using std::string;
using std::vector;

// Enum for mutation types
enum class MutationType {
  SUBSTITUTION, // Base substitution
  INSERTION, // Insertion of bases, to the left of the current position
  DELETION // Deletion of bases
};

// Add operator<< for MutationType
inline std::ostream& operator<<(std::ostream& os, const MutationType& type)
{
  switch (type) {
  case MutationType::SUBSTITUTION:
    os << "SUB";
    break;
  case MutationType::INSERTION:
    os << "INS";
    break;
  case MutationType::DELETION:
    os << "DEL";
    break;
  default:
    os << "UNK";
    break;
  }
  return os;
}

// Structure to represent a mutation
struct Mutation {
  MutationType type; // Type of mutation
  // Position of the mutation on the contig (absolute coordinate)
  uint32_t position;
  string nts; // Bases involved: SUB: src+tgt, INS: inserted, DEL: deleted

  Mutation(MutationType type, uint32_t position, const string& nts = "")
      : type(type)
      , position(position)
      , nts(nts)
  {
  }

  // Create a unique string key for this mutation on a given contig
  std::string create_key(uint32_t contig_index) const;

  string to_string() const
  {
    switch (type) {
    case MutationType::SUBSTITUTION:
      // Ensure nts has at least 2 characters for substitution
      if (nts.length() >= 2) {
        return string(1, nts[0]) + ":" + string(1, nts[1]);
      } else {
        // Handle error or return a default value
        return "ERR_SUB"; // Or throw an exception, log error, etc.
      }
    case MutationType::INSERTION:
      return "+" + nts;
    case MutationType::DELETION:
      return "-" + nts;
    default:
      return "UNK";
    }
  }
};

// Basic data structures for alignment data
struct Contig {
  string id;
  uint32_t length;

  Contig(const string& id = "", uint32_t length = 0)
      : id(id)
      , length(length)
  {
  }
};

struct Read {
  string id;
  uint32_t length;

  Read(const string& id = "", uint32_t length = 0)
      : id(id)
      , length(length)
  {
  }
};

struct Alignment {
  uint32_t read_index;
  uint32_t contig_index;
  uint32_t read_start;
  uint32_t read_end;
  uint32_t contig_start;
  uint32_t contig_end;
  bool is_reverse;
  // Indices into the global mutation store for the corresponding contig
  vector<uint32_t> mutations;

  Alignment(uint32_t read_idx = 0, uint32_t contig_idx = 0,
      uint32_t c_start = 0, uint32_t c_end = 0,
      uint32_t r_start = 0, uint32_t r_end = 0,
      bool is_rev = false)
      : read_index(read_idx)
      , contig_index(contig_idx)
      , read_start(r_start)
      , read_end(r_end)
      , contig_start(c_start)
      , contig_end(c_end)
      , is_reverse(is_rev)
  {
  }

  // Add a mutation index to the alignment
  void add_mutation_index(uint32_t mutation_index)
  {
    mutations.push_back(mutation_index);
  }

  // Clear all mutation indices
  void clear_mutations()
  {
    mutations.clear();
  }
};

// Structure to represent an interval
struct Interval {
  string contig;
  uint32_t start;
  uint32_t end;

  Interval(const string& contig = "", uint32_t start = 0, uint32_t end = 0)
      : contig(contig)
      , start(start)
      , end(end)
  {
  }

  string to_string() const
  {
    return contig + ":" + std::to_string(start) + "-" + std::to_string(end);
  }
};
