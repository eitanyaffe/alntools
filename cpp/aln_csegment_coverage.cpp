#include "CsegmentCovMatrix.h"
#include "Params.h"
#include "utils.h"
#include <iostream>

using namespace std;

void csegment_coverage_params(const char* name, int argc, char** argv, Parameters& params)
{
    params.add_parser("ifn_libraries", new ParserFilename("library table with lib_id and aln_fn columns"), true);
    params.add_parser("ifn_segments", new ParserFilename("segment table with segment, contig, start, end columns"), true);
    params.add_parser("ifn_cseg_map", new ParserFilename("csegment map with segment and csegment columns"), true);
    params.add_parser("ofn_coverage", new ParserFilename("output coverage matrix file with unique read counts per csegment"), true);
    params.add_parser("ofn_read_counts", new ParserFilename("output read counts table per library"), true);
    params.add_parser("min_segment_length", new ParserInteger("minimum segment length for filtering (default 1000)", 1000), false);
    params.add_parser("clip_mode", new ParserString("clipping mode (all, complete, allow_one_side_clip, only_one_side_clipped, only_two_side_clipped, only_clipped, local_align)", "complete"), false);
    params.add_parser("clip_margin", new ParserInteger("clipping margin in bases (default 10)", 10), false);
    params.add_parser("min_mutations_percent", new ParserDouble("minimum mutations percentage (default 0.0)", 0.0), false);
    params.add_parser("max_mutations_percent", new ParserDouble("maximum mutations percentage (default 0.1)", 0.1), false);
    params.add_parser("min_alignment_length", new ParserInteger("minimum alignment length in read coordinates (default 1000)", 1000), false);
    params.add_parser("max_alignment_length", new ParserInteger("maximum alignment length in read coordinates (default 0, no limit)", 0), false);
    params.add_parser("min_indel_length", new ParserInteger("minimum indel length to include in mutation density calculations (default 3)", 3), false);

    if (argc == 1) {
        params.usage(name);
        exit(1);
    }

    params.read(argc, argv);
    params.parse();
    params.verify_mandatory();
    params.print(cout);
}

int csegment_coverage_main(const char* name, int argc, char** argv)
{
    Parameters params;
    csegment_coverage_params(name, argc, argv, params);

    string ifn_libraries = params.get_string("ifn_libraries");
    string ifn_segments = params.get_string("ifn_segments");
    string ifn_cseg_map = params.get_string("ifn_cseg_map");
    string ofn_coverage = params.get_string("ofn_coverage");
    string ofn_read_counts = params.get_string("ofn_read_counts");
    int min_segment_length = params.get_int("min_segment_length");
    
    ClipMode clip_mode = string_to_clip_mode(params.get_string("clip_mode"));
    int clip_margin = params.get_int("clip_margin");
    double min_mutations_percent = params.get_double("min_mutations_percent");
    double max_mutations_percent = params.get_double("max_mutations_percent");
    int min_alignment_length = params.get_int("min_alignment_length");
    int max_alignment_length = params.get_int("max_alignment_length");
    int min_indel_length = params.get_int("min_indel_length");

    CsegmentCovMatrix csegment_cov_matrix;
    csegment_cov_matrix.compute(ifn_libraries, ifn_segments, ifn_cseg_map, ofn_coverage, ofn_read_counts,
                               min_segment_length, clip_mode, clip_margin,
                               min_mutations_percent, max_mutations_percent,
                               min_alignment_length, max_alignment_length, min_indel_length);

    return 0;
}

