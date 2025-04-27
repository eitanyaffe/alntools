#include <cstdlib>
#include <iostream>
#include <string>

#include "Params.h"
#include "alignment_store.h"
#include "paf_reader.h"

using namespace std;

void construct_command(
    const string& ifn_paf,
    const string& ifn_contigs,
    const string& ifn_reads,
    bool should_verify,
    const string& aln_file, int max_reads,
    bool quit_on_error)
{
  PafReader reader;
  AlignmentStore store;

  if (should_verify) {
    cout << "Loading reads and contigs...\n";
    reader.load_reads_contigs(ifn_reads, ifn_contigs);
  }

  cout << "Reading PAF file: " << ifn_paf << "\n";
  reader.read_paf(ifn_paf, store, max_reads, should_verify, quit_on_error);

  cout << "Writing alignment file: " << aln_file << "\n";
  store.save(aln_file);

  cout << "Store info:\n"
       << "  Reads: " << store.get_read_count() << "\n"
       << "  Alignments: " << store.get_alignment_count() << "\n";

  cout << "Done! Processed " << store.get_alignment_count() << " alignments\n";
}

void construct_params(const char* name, int argc, char** argv, Parameters& params)
{
  params.add_parser("ifn_paf", new ParserFilename("input alignment PAF file"), false);
  params.add_parser("ofn", new ParserFilename("output ALN file"), true);
  params.add_parser("verify", new ParserBoolean("should verify PAF file using reads and contigs", false), false);
  params.add_parser("ifn_reads",
      new ParserFilename("input read FASTQ file (used only if verifying alignments)"), false);
  params.add_parser("ifn_contigs",
      new ParserFilename("input contig FASTA file (used only if verifying alignments)"), false);
  params.add_parser("max_reads", new ParserInteger("use only this number of alignments (0: all)", 0), false);
  params.add_parser("quit_on_error", new ParserBoolean("quit on error", true), false);

  if (argc == 1) {
    params.usage(name);
    exit(1);
  }

  // read command line params
  params.read(argc, argv);
  params.parse();
  params.verify_mandatory();
  params.print(cout);
}

int construct_main(const char* name, int argc, char** argv)
{
  Parameters params;
  construct_params(name, argc, argv, params);

  string ifn_paf = params.get_string("ifn_paf");
  string ifn_contigs = params.get_string("ifn_contigs");
  string ifn_reads = params.get_string("ifn_reads");
  bool should_verify = params.get_bool("verify");
  string ofn = params.get_string("ofn");
  int max_reads = params.get_int("max_reads");
  bool quit_on_error = params.get_bool("quit_on_error");
  construct_command(ifn_paf, ifn_contigs, ifn_reads, should_verify, ofn, max_reads, quit_on_error);

  return 0;
}
