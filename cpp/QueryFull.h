#ifndef QUERYFULL_H
#define QUERYFULL_H

#include "alignment_store.h"
#include "aln_types.h"
#include <cstdint>
#include <string>
#include <vector>

struct FullOutputAlignments {
  uint64_t alignment_index;
  std::string read_id;
  std::string contig_id;
  int read_start;
  int read_end;
  int contig_start;
  int contig_end;
  bool is_reverse;
  std::string cs_tag;
};

struct FullOutputMutations {
  uint64_t alignment_index;
  std::string read_id;
  std::string contig_id;
  MutationType type;
  int position;
  std::string desc;
};

class QueryFull {
  private:
  std::vector<Interval> intervals;
  const AlignmentStore& store;

  std::vector<FullOutputAlignments> output_alignments;
  std::vector<FullOutputMutations> output_mutations;

  void generate_output_data();

  public:
  QueryFull(const std::vector<Interval>& intervals,
      const AlignmentStore& store);

  // execute the query
  void execute();

  // write the output rows to a table
  void write_to_csv(const std::string& ofn_prefix);

  // getters
  const std::vector<FullOutputAlignments>& get_output_alignments() const;
  const std::vector<FullOutputMutations>& get_output_mutations() const;
};

#endif // QUERYFULL_H