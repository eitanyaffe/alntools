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

// convert string to HeightStyle enum
HeightStyle string_to_height_style(const std::string& style_str)
{
  if (style_str == "by_mutations") {
    return HeightStyle::BY_MUTATIONS;
  } else {
    // default to BY_COORD
    return HeightStyle::BY_COORD;
  }
}

void query_params(const char* name, int argc, char** argv, Parameters& params)
{
  params.add_parser("ifn_aln", new ParserFilename("input ALN file"), true);
  params.add_parser("ifn_intervals", new ParserFilename("input table with query contig intervals"), true);
  params.add_parser("ofn_prefix", new ParserFilename("output tab-delimited table prefix"), true);
  params.add_parser("mode", new ParserString("query mode (full, pileup, bin)", "full"), true);
  params.add_parser("pileup_mode", new ParserString("pileup report mode (all, covered, mutated)", "covered"), false);
  params.add_parser("binsize", new ParserInteger("bin size for 'bin' mode", 100), false);
  params.add_parser("height_style", new ParserString("alignment height style for 'full' mode (by_coord, by_mutations)", "by_coord"), false);

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

  // Validate height_style if mode is 'full'
  if (mode == "full") {
    string height_style = params.get_string("height_style");
    if (height_style != "by_coord" && height_style != "by_mutations") {
      cerr << "error: invalid height_style specified: " << height_style << ". Must be 'by_coord' or 'by_mutations'." << endl;
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

  // Get pileup mode string and convert to enum
  PileupReportMode pileup_mode = string_to_pileup_report_mode(params.get_string("pileup_mode"));

  // Get height style and convert to enum
  HeightStyle height_style = string_to_height_style(params.get_string("height_style"));

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
  if (mode == "full") {
    cout << "  height_style: " << params.get_string("height_style") << endl;
  }

  vector<Interval> intervals;
  read_intervals(ifn_intervals, intervals);
  cout << "read " << intervals.size() << " intervals from " << ifn_intervals << endl;

  AlignmentStore store;
  store.load(ifn_aln);

  if (mode == "full") {
    QueryFull queryFull(intervals, store, height_style);
    queryFull.execute();
    queryFull.write_to_csv(ofn_prefix);
  } else if (mode == "pileup") {
    QueryPileup queryPileup(intervals, store, pileup_mode);
    queryPileup.execute();
    queryPileup.write_to_csv(ofn_prefix);
  } else if (mode == "bin") {
    QueryBin queryBin(intervals, store, binsize);
    queryBin.execute();
    queryBin.write_to_csv(ofn_prefix);
  }

  return 0;
}