#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "alignment_store.h"
#include "Params.h"
#include "utils.h"
#include <algorithm>

using namespace std;

// Main verification function
void verify_command(const string &ifn_aln,
                    const string &ifn_reads,
                    const string &ifn_contigs,
                    int max_reads,
                    const string &ofn_reads,
                    const string &ofn_contigs)
{
  AlignmentStore store;
  cout << "Reading alignment file: " << ifn_aln << "\n";
  store.load(ifn_aln);

  // limit the number of alignments to verify
  vector<Alignment> alignments = store.get_alignments();
  if (max_reads > 0)
  {
    alignments = vector<Alignment>(alignments.begin(), alignments.begin() + max_reads);
  }

  // Collect contig and read IDs from the store
  vector<string> contig_ids, read_ids;
  for (const auto &alignment : alignments)
  {
    contig_ids.push_back(store.get_contig_id(alignment.contig_index));
    read_ids.push_back(store.get_read_id(alignment.read_index));
  }

  unordered_set<string> contig_set(contig_ids.begin(), contig_ids.end());
  unordered_set<string> read_set(read_ids.begin(), read_ids.end());

  unordered_map<string, string> contigs;
  unordered_map<string, string> reads;
  read_fasta(ifn_contigs, contig_set, contigs);
  read_fastq(ifn_reads, read_set, reads);

  write_fasta(ofn_contigs, contigs);
  write_fastq(ofn_reads, reads);

  int bad_alignment_count = 0;
  for (const auto &alignment : alignments)
  {

    const string &contig_id = store.get_contig_id(alignment.contig_index);
    const string &read_id = store.get_read_id(alignment.read_index);

    cout << "==================\n"
         << "Read: " << read_id << " [" << alignment.read_start << "," << alignment.read_end << "] "
         << "  Contig: " << contig_id << " [" << alignment.contig_start << "," << alignment.contig_end << "] "
         << "  Is reverse: " << (alignment.is_reverse ? "yes" : "no") << "\n";

    // Count mutations by type
    auto count_mutations_by_type = [](const vector<Mutation> &mutations)
    {
      size_t subs = 0, ins = 0, dels = 0;
      for (const auto &mut : mutations)
      {
        switch (mut.type)
        {
        case MutationType::SUBSTITUTION:
          subs++;
          break;
        case MutationType::INSERTION:
          ins++;
          break;
        case MutationType::DELETION:
          dels++;
          break;
        }
      }
      return make_tuple(subs, ins, dels);
    };

    auto [num_subs, num_ins, num_dels] = count_mutations_by_type(alignment.mutations);
    cout << "Mutations - Substitutions: " << num_subs
         << ", Insertions: " << num_ins
         << ", Deletions: " << num_dels << "\n";

    massert(contigs.find(contig_id) != contigs.end(),
            "Error: Contig '%s' not found in FASTA file.", contig_id.c_str());
    massert(reads.find(read_id) != reads.end(),
            "Error: Read '%s' not found in FASTQ file.", read_id.c_str());

    cout << "mutating contig with " << alignment.mutations.size() << " mutations" << endl;

    massert(contigs.find(contig_id) != contigs.end(), "contig %s not found in map", contig_id.c_str());
    massert(reads.find(read_id) != reads.end(), "read %s not found in map", read_id.c_str());

    string contig_fragment = contigs[contig_id].substr(alignment.contig_start,
                                                       alignment.contig_end - alignment.contig_start);

    string mutated_contig = apply_mutations(contig_fragment, alignment.mutations);
    string read_segment = reads[read_id].substr(alignment.read_start,
                                                alignment.read_end - alignment.read_start);

    if (alignment.is_reverse)
      mutated_contig = reverse_complement(mutated_contig);

    massert(read_segment.size() == mutated_contig.size(),
            "read segment length (%zu) does not match mutated contig length (%zu)",
            read_segment.size(), mutated_contig.size());

    bool mismatch_found = false;
    for (size_t i = 0; i < read_segment.size(); ++i)
    {
      if (mutated_contig[i] != read_segment[i])
      {

        // Calculate start and end indices for the segment
        size_t start = (i >= 5) ? i - 8 : 0;
        size_t end = (i + 5 < read_segment.size()) ? i + 8 : read_segment.size() - 1;

        cout << "Mismatch found, fragment coordinate=" << i << endl;
        cout << "read        : " << read_segment.substr(start, end - start + 1) << "\n";
        cout << "contig_mut  : " << mutated_contig.substr(start, end - start + 1) << "\n";
        cout << "contig_orig : " << contig_fragment.substr(start, end - start + 1) << "\n";

        mismatch_found = true;
        break;
      }
    }

    if (mismatch_found)
    {
      bad_alignment_count++;
      if (bad_alignment_count > 100)
      {
        cerr << "Too many bad alignments. Exiting.\n";
        exit(-1);
      }
    }
    else
    {
      cout << "Alignment is good.\n";
    }
  }

  cout << "Verification complete. Total alignments processed: " << store.get_alignment_count() << "\n";
  cout << "Bad alignments found: " << bad_alignment_count << " out of " << store.get_alignment_count()
       << " (" << (bad_alignment_count * 100.0 / store.get_alignment_count()) << "%)\n";
}

void verify_params(const char *name, int argc, char **argv, Parameters &params)
{
  params.add_parser("ifn_aln", new ParserFilename("input ALN file"), true);
  params.add_parser("ifn_reads", new ParserFilename("input reads, FASTQ"), true);
  params.add_parser("ifn_contigs", new ParserFilename("input contigs, FASTA"), true);

  params.add_parser("max_reads", new ParserInteger("use only this number of alignments (0: all)", 0), false);

  params.add_parser("ofn_contigs", new ParserFilename("contigs limited to alignments, FASTA"), false);
  params.add_parser("ofn_reads", new ParserFilename("reads limited to alignments, FASTA"), false);

  if (argc == 1)
  {
    params.usage(name);
    exit(1);
  }

  // read command line params
  params.read(argc, argv);
  params.parse();
  params.verify_mandatory();
  params.print(cout);
}

int verify_main(const char *name, int argc, char **argv)
{
  Parameters params;
  verify_params(name, argc, argv, params);

  string ifn_aln = params.get_string("ifn_aln");
  string ifn_reads = params.get_string("ifn_reads");
  string ifn_contigs = params.get_string("ifn_contigs");

  int max_reads = params.get_int("max_reads");

  string ofn_reads = params.get_string("ofn_reads");
  string ofn_contigs = params.get_string("ofn_contigs");

  verify_command(ifn_aln, ifn_reads, ifn_contigs, max_reads, ofn_reads, ofn_contigs);

  return 0;
}
