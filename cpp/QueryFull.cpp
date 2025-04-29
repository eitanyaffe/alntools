#include "QueryFull.h"
#include "utils.h"
#include <fstream>
#include <iostream>

using namespace std;

QueryFull::QueryFull(const vector<Interval>& intervals, const string& ofn_prefix, const AlignmentStore& store)
    : intervals(intervals)
    , ofn_prefix(ofn_prefix)
    , store(store)
{
}

void QueryFull::write_to_csv()
{
  cout << "writing to " << ofn_prefix + "_alignments.tsv" << endl;
  ofstream ofs(ofn_prefix + "_alignments.tsv");

  // check if file opened successfully
  if (!ofs.is_open()) {
    cerr << "error: could not open file " << ofn_prefix + "_alignments.tsv" << endl;
    exit(1);
  }

  int alignment_count = 0;
  // Write header
  ofs << "read_id\tcontig_id\tread_start\tread_end\tcontig_start\tcontig_end\tis_reverse\tcs_tag\n";

  for (const auto& interval : intervals) {
    std::vector<std::reference_wrapper<const Alignment>> alignments = store.get_alignments_in_interval(interval);
    for (const auto& alignment : alignments) {
      const auto& aln = alignment.get();
      string read_id = store.get_read_id(aln.read_index);
      string contig_id = store.get_contig_id(aln.contig_index);
      string cs_string = generate_cs_tag(aln);
      ofs << read_id << "\t"
          << contig_id << "\t"
          << aln.read_start << "\t"
          << aln.read_end << "\t"
          << aln.contig_start << "\t"
          << aln.contig_end << "\t"
          << (aln.is_reverse ? "true" : "false") << "\t"
          << cs_string << "\n";
      alignment_count++;
    }
  }
  ofs.close();
  cout << "wrote " << alignment_count << " alignments to " << ofn_prefix + "_alignments.tsv" << endl;

  // output mutations
  cout << "writing to " << ofn_prefix + "_mutations.tsv" << endl;
  ofstream ofs_mutations(ofn_prefix + "_mutations.tsv");
  // check if file opened successfully
  if (!ofs_mutations.is_open()) {
    cerr << "error: could not open file " << ofn_prefix + "_mutations.tsv" << endl;
    exit(1);
  }

  int mutation_count = 0;
  ofs_mutations << "read_id\tcontig_id\tread_start\tread_end\tcontig_start\tcontig_end\tis_reverse\tmutation_type\tmutation_position\tmutation_desc\n";
  for (const auto& interval : intervals) {
    std::vector<std::reference_wrapper<const Alignment>> alignments = store.get_alignments_in_interval(interval);
    for (const auto& alignment : alignments) {
      const auto& aln = alignment.get();
      string read_id = store.get_read_id(aln.read_index);
      string contig_id = store.get_contig_id(aln.contig_index);
      for (const auto& mutation : aln.mutations) {
        ofs_mutations << read_id << "\t"
                      << contig_id << "\t"
                      << mutation.type << "\t"
                      << mutation.position << "\t"
                      << mutation.to_string() << "\n";
        mutation_count++;
      }
    }
  }
  ofs_mutations.close();
  cout << "wrote " << mutation_count << " mutations to " << ofn_prefix + "_mutations.tsv" << endl;
}
