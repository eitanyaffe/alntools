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
  uint32_t position; // Position relative to the start of the alignment
  string read_nts; // Bases in the query sequence (for substitutions and insertions)
  string ref_nts; // Bases in the target sequence (for substitutions and deletions)

  Mutation(MutationType type, uint32_t position, const string& read_nts = "",
      const string& ref_nts = "")
      : type(type)
      , position(position)
      , read_nts(read_nts)
      , ref_nts(ref_nts)
  {
  }

  string to_string() const
  {
    switch (type) {
    case MutationType::SUBSTITUTION:
      return read_nts + ":" + ref_nts;
    case MutationType::INSERTION:
      return "+" + read_nts;
    case MutationType::DELETION:
      return "-" + ref_nts;
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
  vector<Mutation> mutations; // Vector to store mutations

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

  // Add a mutation to the alignment
  void add_mutation(const Mutation& mutation)
  {
    mutations.push_back(mutation);
  }

  // Clear all mutations
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
};
