#include <iostream>
#include <string>
#include <vector>
#include "Params.h" // Assuming Params class is defined here
#include "utils.h"	// Assuming massert is defined here
#include "alignment_store.h"
using namespace std;

void query_full_command(const vector<Interval> &intervals,
						const string &ofn, const AlignmentStore &store)
{
	cout << "query full command called" << endl;
	ofstream ofs(ofn);
	
	// Write header
	ofs << "read_id\tcontig_id\tread_start\tread_end\tcontig_start\tcontig_end\tis_reverse\tcs_tag\n";
	
	for (const auto &interval : intervals)
	{
		std::vector<std::reference_wrapper<const Alignment>> alignments = store.get_alignments_in_interval(interval);
		for (const auto &alignment : alignments)
		{
			const auto &aln = alignment.get();
			string read_id = store.get_read_id(aln.read_index);
			string contig_id = store.get_contig_id(aln.contig_index);
			
			ofs << read_id << "\t"
				<< contig_id << "\t" 
				<< aln.read_start << "\t"
				<< aln.read_end << "\t"
				<< aln.contig_start << "\t"
				<< aln.contig_end << "\t"
				<< (aln.is_reverse ? "true" : "false") << "\t";
			ofs << generate_cs_tag(aln) << "\n";
		}
	}
}

void query_pileup_command(const vector<Interval> &intervals,
						  const string &ofn, const AlignmentStore &store)
{
	cout << "query pileup command called" << endl;
}

void query_bin_command(const vector<Interval> &intervals,
					   const string &ofn, const AlignmentStore &store, int binsize)
{
	cout << "query bin command called" << endl;
}

void query_params(const char *name, int argc, char **argv, Parameters &params)
{
	params.add_parser("ifn_aln", new ParserFilename("input ALN file"), true);
	params.add_parser("ifn_intervals", new ParserFilename("input table with query contig intervals"), true);
	params.add_parser("ofn", new ParserFilename("output tab-delimited table"), true);
	params.add_parser("mode", new ParserString("query mode (full, pileup, bin)", "full"), true); // Defaulting to "full" for now, update if needed
	params.add_parser("binsize", new ParserInteger("bin size for 'bin' mode", 100), false);

	if (argc == 1)
	{
		params.usage(name);
		exit(1);
	}

	// read command line params
	params.read(argc, argv);
	params.parse();
	params.verify_mandatory();

	// Validate mode
	string mode = params.get_string("mode");
	if (mode != "full" && mode != "pileup" && mode != "bin")
	{
		cerr << "error: invalid mode specified: " << mode << ". Must be 'full', 'pileup', or 'bin'." << endl;
		exit(1);
	}

	// If mode is 'bin', binsize must be positive
	if (mode == "bin")
	{
		int binsize = params.get_int("binsize");
		if (binsize <= 0)
		{
			cerr << "error: binsize must be a positive integer for mode 'bin'." << endl;
			exit(1);
		}
	}

	params.print(cout);
}

int query_main(const char *name, int argc, char **argv)
{
	Parameters params;
	query_params(name, argc, argv, params);

	string ifn_aln = params.get_string("ifn_aln");
	string ifn_intervals = params.get_string("ifn_intervals");
	string ofn = params.get_string("ofn");
	string mode = params.get_string("mode");
	int binsize = params.get_int("binsize"); // Will be 0 if not specified or mode is not 'bin'

	cout << "query command called:" << endl;
	cout << "  ifn_aln: " << ifn_aln << endl;
	cout << "  ifn_intervals: " << ifn_intervals << endl;
	cout << "  ofn: " << ofn << endl;
	cout << "  mode: " << mode << endl;
	if (mode == "bin")
	{
		cout << "  binsize: " << binsize << endl;
	}

	vector<Interval> intervals;
	read_intervals(ifn_intervals, intervals);
	cout << "read " << intervals.size() << " intervals from " << ifn_intervals << endl;

	AlignmentStore store;
	store.load(ifn_aln);

	if (mode == "full")
	{
		query_full_command(intervals, ofn, store);
	}
	else if (mode == "pileup")
	{
		query_pileup_command(intervals, ofn, store);
	}
	else if (mode == "bin")
	{
		query_bin_command(intervals, ofn, store, binsize);
	}

	return 0;
}