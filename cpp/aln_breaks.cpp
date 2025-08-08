#include <fstream>
#include <iostream>
#include <string>

#include "Params.h"
#include "alignment_store.h"
#include "utils.h"

using namespace std;

void breaks_command(const string& aln_file, const string& output_file, 
                   uint32_t window_size, double p_threshold, uint32_t min_reads)
{
  cout << "loading alignment file " << aln_file << endl;

  AlignmentStore store;
  store.load(aln_file);

  cout << "finding break positions with window size " << window_size 
       << ", p threshold " << p_threshold 
       << ", and min reads " << min_reads << endl;

  auto results = store.find_break_positions(window_size, p_threshold, min_reads);

  cout << "writing results to " << output_file << endl;
  ofstream out(output_file);
  massert(out.is_open(), "failed to open file for writing: %s", output_file.c_str());

  // write header
  out << "contig\tposition\torientation\tt\te\tenrichment\tpval\tqval\n";

  // write results
  for (const auto& result : results) {
    out << result.contig_id << "\t"
        << result.position << "\t"
        << result.orientation << "\t"
        << result.t << "\t"
        << result.e << "\t"
        << result.enrichment << "\t"
        << result.pval << "\t"
        << result.qval << "\n";
  }

  out.close();

  cout << "found " << results.size() << " significant break positions" << endl;
}

void breaks_params(const char* name, int argc, char** argv, Parameters& params)
{
  params.add_parser("ifn", new ParserFilename("input ALN file"), true);
  params.add_parser("ofn", new ParserFilename("output TSV file"), true);
  params.add_parser("window", new ParserInteger("window size for background calculation"), true);
  params.add_parser("pval", new ParserDouble("p-value threshold for reporting"), true);
  params.add_parser("min_reads", new ParserInteger("minimum reads required to test position"), false);

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

int breaks_main(const char* name, int argc, char** argv)
{
  Parameters params;
  breaks_params(name, argc, argv, params);

  string ifn = params.get_string("ifn");
  string ofn = params.get_string("ofn");
  uint32_t window_size = params.get_int("window");
  double p_threshold = params.get_double("pval");
  uint32_t min_reads = params.is_defined("min_reads") ? params.get_int("min_reads") : 1;

  breaks_command(ifn, ofn, window_size, p_threshold, min_reads);

  return 0;
}