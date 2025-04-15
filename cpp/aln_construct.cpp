#include <iostream>
#include <string>
#include <cstdlib>

#include "Params.h"
#include "paf_reader.h"
#include "alignment_store.h"

using namespace std;

void construct_command(const string& paf_file, const string& aln_file) 
{
    cout << "Reading PAF file: " << paf_file << "\n";
    
    PafReader reader;
    AlignmentStore store;
    
    cout << "Processing alignments...\n";
    reader.read_paf(paf_file, store);
    
    cout << "Writing alignment file: " << aln_file << "\n";
    store.save(aln_file);
    
    cout << "Store info:\n"
         << "  Reads: " << store.get_read_count() << "\n"
         << "  Alignments: " << store.get_alignment_count() << "\n";

    cout << "Done! Processed " << store.get_alignment_count() << " alignments\n";
}

void construct_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn", new ParserFilename("input PAF file"), true);
  params.add_parser("ofn", new ParserFilename("output ALN file"), true);
  
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

int construct_main(const char* name, int argc, char **argv)
{
  Parameters params;
  construct_params(name, argc, argv, params);
  
  string ifn = params.get_string("ifn");
  string ofn = params.get_string("ofn");

  construct_command(ifn, ofn);
  
  return 0;
}
