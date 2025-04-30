#include "QueryFull.h"
#include "utils.h"
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

QueryFull::QueryFull(const vector<Interval>& intervals, const AlignmentStore& store)
    : intervals(intervals)
    , store(store)
{
}

void QueryFull::generate_output_data()
{
  uint64_t current_alignment_index = 0;
  output_alignments.clear();
  output_mutations.clear();

  for (const auto& interval : intervals) {
    std::vector<std::reference_wrapper<const Alignment>> alignments = store.get_alignments_in_interval(interval);
    for (const auto& alignment_ref : alignments) {
      const auto& aln = alignment_ref.get();
      string read_id = store.get_read_id(aln.read_index);
      string contig_id = store.get_contig_id(aln.contig_index);
      string cs_string = generate_cs_tag(aln);

      output_alignments.push_back({ current_alignment_index,
          read_id,
          contig_id,
          static_cast<int>(aln.read_start),
          static_cast<int>(aln.read_end),
          static_cast<int>(aln.contig_start),
          static_cast<int>(aln.contig_end),
          aln.is_reverse,
          cs_string });

      for (const auto& mutation : aln.mutations) {
        output_mutations.push_back({ current_alignment_index,
            read_id,
            contig_id,
            mutation.type,
            static_cast<int>(mutation.position),
            mutation.to_string() });
      }
      current_alignment_index++;
    }
  }
}

void QueryFull::write_to_csv(const std::string& ofn_prefix)
{
  // --- Write Alignments ---
  cout << "writing alignments to " << ofn_prefix + "_alignments.tsv" << endl;
  ofstream ofs_alignments(ofn_prefix + "_alignments.tsv");
  if (!ofs_alignments.is_open()) {
    cerr << "error: could not open file " << ofn_prefix + "_alignments.tsv" << endl;
    exit(1);
  }

  // Write header
  ofs_alignments << "alignment_index\tread_id\tcontig_id\tread_start\tread_end\tcontig_start\tcontig_end\tis_reverse\tcs_tag\n"; // Corrected header

  // Write data
  for (const auto& aln_data : output_alignments) {
    ofs_alignments << aln_data.alignment_index << "\t"
                   << aln_data.read_id << "\t"
                   << aln_data.contig_id << "\t"
                   << aln_data.read_start << "\t"
                   << aln_data.read_end << "\t"
                   << aln_data.contig_start << "\t"
                   << aln_data.contig_end << "\t"
                   << (aln_data.is_reverse ? "true" : "false") << "\t"
                   << aln_data.cs_tag << "\n"; // Corrected data line format
  } // Moved brace to its own line
  ofs_alignments.close();
  cout << "wrote " << output_alignments.size() << " alignments to " << ofn_prefix + "_alignments.tsv" << endl;

  // --- Write Mutations ---
  cout << "writing mutations to " << ofn_prefix + "_mutations.tsv" << endl;
  ofstream ofs_mutations(ofn_prefix + "_mutations.tsv");
  if (!ofs_mutations.is_open()) {
    cerr << "error: could not open file " << ofn_prefix + "_mutations.tsv" << endl;
    exit(1); // Consider using exceptions or error codes instead of exit(1)
  }

  // Write header
  ofs_mutations << "alignment_index\tread_id\tcontig_id\tmutation_type\tmutation_position\tmutation_desc\n"; // Corrected header

  // Write data
  for (const auto& mut_data : output_mutations) {
    // Assuming MutationType can be streamed directly or needs conversion
    // If conversion needed, e.g., mutation_type_to_string(mut_data.type)
    ofs_mutations << mut_data.alignment_index << "\t"
                  << mut_data.read_id << "\t"
                  << mut_data.contig_id << "\t"
                  << mut_data.type << "\t" // Adjust if type needs string conversion
                  << mut_data.position << "\t"
                  << mut_data.desc << "\n"; // Corrected data line format
  }
  ofs_mutations.close();
  cout << "wrote " << output_mutations.size() << " mutations to " << ofn_prefix + "_mutations.tsv" << endl;
}

void QueryFull::execute()
{
  generate_output_data();
}

const std::vector<FullOutputAlignments>& QueryFull::get_output_alignments() const
{
  return output_alignments;
}

const std::vector<FullOutputMutations>& QueryFull::get_output_mutations() const
{
  return output_mutations;
}
