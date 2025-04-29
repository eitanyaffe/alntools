#include "Params.h"
#include "QueryBin.h"
#include "QueryFull.h"
#include "QueryPileup.h"
#include "alignment_store.h"
#include "utils.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;

void query_params(const char* name, int argc, char** argv, Parameters& params)
{
  params.add_parser("ifn_aln", new ParserFilename("input ALN file"), true);
  params.add_parser("ifn_intervals", new ParserFilename("input table with query contig intervals"), true);
  params.add_parser("ofn_prefix", new ParserFilename("output tab-delimited table prefix"), true);
  params.add_parser("mode", new ParserString("query mode (full, pileup, bin)", "full"), true);
  params.add_parser("pileup_mode", new ParserString("pileup report mode (all, covered, mutated)", "covered"), false);
  params.add_parser("skip_empty_bins", new ParserBoolean("skip bins with no mutations", false), false);
  params.add_parser("binsize", new ParserInteger("bin size for 'bin' mode", 100), false);

  if (argc == 1) {
    params.usage(name);
    exit(1);
  }

  // read command line params
  params.read(argc, argv);
  params.parse();
  params.verify_mandatory();

  // Validate mode
  string mode = params.get_string("mode");
  if (mode != "full" && mode != "pileup" && mode != "bin") {
    cerr << "error: invalid mode specified: " << mode << ". Must be 'full', 'pileup', or 'bin'." << endl;
    exit(1);
  }

  // If mode is 'bin', binsize must be positive
  if (mode == "bin") {
    int binsize = params.get_int("binsize");
    if (binsize <= 0) {
      cerr << "error: binsize must be a positive integer for mode 'bin'." << endl;
      exit(1);
    }
  }

  // If mode is 'pileup', validate pileup_mode
  if (mode == "pileup") {
    string pileup_mode_str = params.get_string("pileup_mode");
    if (pileup_mode_str != "all" && pileup_mode_str != "covered" && pileup_mode_str != "mutated") {
      cerr << "error: invalid pileup_mode specified: " << pileup_mode_str
           << ". Must be 'all', 'covered', or 'mutated'." << endl;
      exit(1);
    }
  }

  params.print(cout);
}

int query_main(const char* name, int argc, char** argv)
{
  Parameters params;
  query_params(name, argc, argv, params);

  string ifn_aln = params.get_string("ifn_aln");
  string ifn_intervals = params.get_string("ifn_intervals");
  string ofn_prefix = params.get_string("ofn_prefix");
  string mode = params.get_string("mode");
  int binsize = params.get_int("binsize"); // Will be 0 if not specified or mode is not 'bin'
  bool skip_empty_bins = params.get_bool("skip_empty_bins");

  // Get pileup mode string and convert to enum
  PileupReportMode pileup_mode = PileupReportMode::COVERED; // Default
  if (mode == "pileup") {
    string pileup_mode_str = params.get_string("pileup_mode");
    if (pileup_mode_str == "all") {
      pileup_mode = PileupReportMode::ALL;
    } else if (pileup_mode_str == "mutated") {
      pileup_mode = PileupReportMode::MUTATED;
    } // else it remains COVERED (the default)
  }

  cout << "query command called:" << endl;
  cout << "  ifn_aln: " << ifn_aln << endl;
  cout << "  ifn_intervals: " << ifn_intervals << endl;
  cout << "  ofn_prefix: " << ofn_prefix << endl;
  cout << "  mode: " << mode << endl;
  if (mode == "bin") {
    cout << "  binsize: " << binsize << endl;
  }
  if (mode == "pileup") {
    cout << "  pileup_mode: " << params.get_string("pileup_mode") << endl;
  }

  vector<Interval> intervals;
  read_intervals(ifn_intervals, intervals);
  cout << "read " << intervals.size() << " intervals from " << ifn_intervals << endl;

  AlignmentStore store;
  store.load(ifn_aln);

  if (mode == "full") {
    QueryFull queryFull(intervals, ofn_prefix, store);
    queryFull.write_to_csv();
  } else if (mode == "pileup") {
    QueryPileup queryPileup(intervals, ofn_prefix, store, pileup_mode);
    queryPileup.write_to_csv();
  } else if (mode == "bin") {
    QueryBin queryBin(intervals, ofn_prefix, store, binsize, skip_empty_bins);
    queryBin.write_to_csv();
  }

  return 0;
}