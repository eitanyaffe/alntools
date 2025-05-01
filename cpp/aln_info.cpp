#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>

#include "Params.h"
#include "alignment_store.h"

using namespace std;

void info_command(const string& aln_file)
{
  double size_mb = get_file_size_mb(aln_file);
  cout << "loading alignment file " << aln_file << " (" << size_mb << " MB)" << endl;

  AlignmentStore store;
  store.load(aln_file);

  cout << "Total alignments: " << store.get_alignment_count() << "\n";
  cout << "Total reads: " << store.get_read_count() << "\n";

  // Calculate average alignment length
  size_t total_length = 0;
  for (const auto& aln : store.get_alignments()) {
    total_length += (aln.read_end - aln.read_start);
  }
  double avg_length = store.get_alignment_count() > 0 ? static_cast<double>(total_length) / store.get_alignment_count() : 0;

  cout << "Average alignment length: " << avg_length << " bp\n";

  // Calculate mutation statistics
  size_t total_mutations = 0;
  for (const auto& aln : store.get_alignments()) {
    total_mutations += aln.mutations.size();
  }

  double avg_mutations = store.get_alignment_count() > 0 ? static_cast<double>(total_mutations) / store.get_alignment_count() : 0;

  cout << "Total mutations: " << total_mutations << "\n";
  cout << "Average mutations per alignment: " << avg_mutations << "\n";
}

void info_params(const char* name, int argc, char** argv, Parameters& params)
{
  params.add_parser("ifn", new ParserFilename("input PAF file"), true);

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

int info_main(const char* name, int argc, char** argv)
{
  Parameters params;
  info_params(name, argc, argv, params);

  string ifn = params.get_string("ifn");

  info_command(ifn);

  return 0;
}
