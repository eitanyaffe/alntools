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
  } else if (style_str == "by_coord_left") {
    return HeightStyle::BY_COORD_LEFT;
  } else if (style_str == "by_coord_right") {
    return HeightStyle::BY_COORD_RIGHT;
  } else {
    // default to BY_COORD_LEFT
    return HeightStyle::BY_COORD_LEFT;
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
  params.add_parser("seg_threshold", new ParserDouble("segregating sites threshold for 'bin' mode", 0.2), false);
  params.add_parser("non_ref_threshold", new ParserDouble("non-reference sites threshold for 'bin' mode", 0.9), false);
  params.add_parser("height_style", new ParserString("read height style for 'full' mode (by_coord_left, by_coord_right, by_mutations)", "by_coord_left"), false);
  params.add_parser("max_alignments", new ParserInteger("maximum number of alignments to return per interval (0 for no limit)", 0), false);
  params.add_parser("alignment_filter", new ParserString("alignment filter for 'full' mode (all, single, single_complete, only_multiple)", "all"), false);
  params.add_parser("min_mutations_density", new ParserDouble("minimum mutations density (mutations per 1000bp) for 'full' mode (0.0 for no filter)", 0.0), false);
  params.add_parser("num_threads", new ParserInteger("number of threads for 'bin' mode (0 for auto)", 0), false);

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
    if (height_style != "by_coord_left" && height_style != "by_coord_right" && 
        height_style != "by_mutations") {
      cerr << "error: invalid height_style specified: " << height_style << ". Must be 'by_coord_left', 'by_coord_right', or 'by_mutations'." << endl;
      exit(1);
    }
    
    // validate alignment_filter if mode is 'full'
    string alignment_filter = params.get_string("alignment_filter");
    if (alignment_filter != "all" && alignment_filter != "single" && 
        alignment_filter != "single_complete" && alignment_filter != "only_multiple") {
      cerr << "error: invalid alignment_filter specified: " << alignment_filter << ". Must be 'all', 'single', 'single_complete', or 'only_multiple'." << endl;
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
  double seg_threshold = params.get_double("seg_threshold");
  double non_ref_threshold = params.get_double("non_ref_threshold");

  // Get pileup mode string and convert to enum
  PileupReportMode pileup_mode = string_to_pileup_report_mode(params.get_string("pileup_mode"));

  // Get height style and convert to enum
  HeightStyle height_style = string_to_height_style(params.get_string("height_style"));
  
  // get max alignments parameter
  int max_alignments = params.get_int("max_alignments");
  
  // get alignment filter parameter
  string alignment_filter = params.get_string("alignment_filter");
  
  // get min mutations density parameter
  double min_mutations_density = params.get_double("min_mutations_density");

  cout << "query command called:" << endl;
  cout << "  ifn_aln: " << ifn_aln << endl;
  cout << "  ifn_intervals: " << ifn_intervals << endl;
  cout << "  ofn_prefix: " << ofn_prefix << endl;
  cout << "  mode: " << mode << endl;
  if (mode == "bin") {
    cout << "  binsize: " << binsize << endl;
    cout << "  seg_threshold: " << seg_threshold << endl;
    cout << "  non_ref_threshold: " << non_ref_threshold << endl;
  }
  if (mode == "pileup") {
    cout << "  pileup_mode: " << params.get_string("pileup_mode") << endl;
  }
  if (mode == "full") {
    cout << "  height_style: " << params.get_string("height_style") << endl;
    cout << "  max_alignments: " << max_alignments << endl;
    cout << "  alignment_filter: " << alignment_filter << endl;
    cout << "  min_mutations_density: " << min_mutations_density << endl;
  }

  vector<Interval> intervals;
  read_intervals(ifn_intervals, intervals);
  cout << "read " << intervals.size() << " intervals from " << ifn_intervals << endl;

  AlignmentStore store;
  store.load(ifn_aln);

  if (mode == "full") {
    QueryFull queryFull(intervals, store, height_style, max_alignments, alignment_filter, min_mutations_density);
    queryFull.execute();
    queryFull.write_to_csv(ofn_prefix);
  } else if (mode == "pileup") {
    QueryPileup queryPileup(intervals, store, pileup_mode);
    queryPileup.execute();
    queryPileup.write_to_csv(ofn_prefix);
  } else if (mode == "bin") {
    int num_threads = params.get_int("num_threads");
    QueryBin queryBin(intervals, store, binsize, seg_threshold, non_ref_threshold, num_threads);
    queryBin.execute();
    queryBin.write_to_csv(ofn_prefix);
  }

  return 0;
}