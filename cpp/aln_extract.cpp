#include <cstdlib>
#include <iostream>
#include <string>

#include "Params.h"
#include "alignment_store.h"

using namespace std;

void extract_command(const string& aln_file, const string& output_prefix)
{
  cout << "Reading alignment file: " << aln_file << "\n";

  AlignmentStore store;
  store.load(aln_file);

  store.export_tab_delimited(output_prefix);
}

void extract_params(const char* name, int argc, char** argv, Parameters& params)
{
  params.add_parser("ifn", new ParserFilename("input ALN file"), true);
  params.add_parser("ofn_prefix", new ParserFilename("output prefix for tables"), true);

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

int extract_main(const char* name, int argc, char** argv)
{
  Parameters params;
  extract_params(name, argc, argv, params);

  string ifn = params.get_string("ifn");
  string ofn_prefix = params.get_string("ofn_prefix");

  extract_command(ifn, ofn_prefix);

  return 0;
}
