#include "GetReadIds.h"
#include "Params.h"
#include "utils.h"
#include <iostream>

using namespace std;

void get_read_ids_params(const char* name, int argc, char** argv, Parameters& params)
{
    params.add_parser("ifn_aln",      new ParserFilename("input ALN file"), true);
    params.add_parser("ifn_segments", new ParserFilename("segment table with contig, start, end, bin_id columns"), true);
    params.add_parser("ofn",          new ParserFilename("output tab-delimited file with read_id and bin_id columns"), true);
    params.add_parser("clip_mode",    new ParserString("clipping mode (all, complete, allow_one_side_clip, only_one_side_clipped, only_two_side_clipped, only_clipped, local_align)", "complete"), false);
    params.add_parser("clip_margin",  new ParserInteger("clipping margin in bases (default 10)", 10), false);
    params.add_parser("min_mutations_percent", new ParserDouble("minimum mutations percentage (default 0.0)", 0.0), false);
    params.add_parser("max_mutations_percent", new ParserDouble("maximum mutations percentage (default 0.1)", 0.1), false);
    params.add_parser("min_alignment_length",  new ParserInteger("minimum alignment length in read coordinates (default 1000)", 1000), false);
    params.add_parser("max_alignment_length",  new ParserInteger("maximum alignment length in read coordinates (default 0, no limit)", 0), false);
    params.add_parser("min_indel_length",      new ParserInteger("minimum indel length to include in mutation density calculations (default 3)", 3), false);

    if (argc == 1) {
        params.usage(name);
        exit(1);
    }

    params.read(argc, argv);
    params.parse();
    params.verify_mandatory();
    params.print(cout);
}

int get_read_ids_main(const char* name, int argc, char** argv)
{
    Parameters params;
    get_read_ids_params(name, argc, argv, params);

    string ifn_aln      = params.get_string("ifn_aln");
    string ifn_segments = params.get_string("ifn_segments");
    string ofn          = params.get_string("ofn");

    ClipMode clip_mode             = string_to_clip_mode(params.get_string("clip_mode"));
    int      clip_margin           = params.get_int("clip_margin");
    double   min_mutations_percent = params.get_double("min_mutations_percent");
    double   max_mutations_percent = params.get_double("max_mutations_percent");
    int      min_alignment_length  = params.get_int("min_alignment_length");
    int      max_alignment_length  = params.get_int("max_alignment_length");
    int      min_indel_length      = params.get_int("min_indel_length");

    GetReadIds get_read_ids;
    get_read_ids.compute(ifn_aln, ifn_segments, ofn,
                         clip_mode, clip_margin,
                         min_mutations_percent, max_mutations_percent,
                         min_alignment_length, max_alignment_length,
                         min_indel_length);

    return 0;
}
