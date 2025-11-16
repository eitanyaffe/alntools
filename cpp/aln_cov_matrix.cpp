#include "CovMatrix.h"
#include "Params.h"
#include "utils.h"
#include <iostream>

using namespace std;

void cov_matrix_params(const char* name, int argc, char** argv, Parameters& params)
{
    params.add_parser("ifn_libraries", new ParserFilename("library table with lib_id and aln_fn columns"), true);
    params.add_parser("ifn_segments", new ParserFilename("segment table with segment, contig, start, end columns"), true);
    params.add_parser("ifn_fasta", new ParserFilename("contig fasta file"), true);
    params.add_parser("ofn_mat", new ParserFilename("output coverage matrix file"), true);
    params.add_parser("ofn_fasta", new ParserFilename("output segment fasta file"), true);
    params.add_parser("actual_nts", new ParserBoolean("use actual nucleotides when writing fasta (default T)", true), false);
    params.add_parser("should_create_fasta", new ParserBoolean("write fasta output (default T)", true), false);
    params.add_parser("min_segment_length", new ParserInteger("minimum segment length for filtering (default 1000)", 1000), false);
    params.add_parser("clip_mode", new ParserString("clipping mode (all, complete, allow_one_side_clip, only_one_side_clipped, only_two_side_clipped, only_clipped, local_align)", "complete"), false);
    params.add_parser("clip_margin", new ParserInteger("clipping margin in bases (default 10)", 10), false);
    params.add_parser("min_mutations_percent", new ParserDouble("minimum mutations percentage (default 0.0)", 0.0), false);
    params.add_parser("max_mutations_percent", new ParserDouble("maximum mutations percentage (default 0.1)", 0.1), false);
    params.add_parser("min_alignment_length", new ParserInteger("minimum alignment length in read coordinates (default 1000)", 1000), false);
    params.add_parser("max_alignment_length", new ParserInteger("maximum alignment length in read coordinates (default 0, no limit)", 0), false);
    params.add_parser("min_indel_length", new ParserInteger("minimum indel length to include in mutation density calculations (default 3)", 3), false);
    params.add_parser("ofn_lib_map", new ParserFilename("output library index to library ID mapping file"), false);

    if (argc == 1) {
        params.usage(name);
        exit(1);
    }

    params.read(argc, argv);
    params.parse();
    params.verify_mandatory();
    params.print(cout);
}

int cov_matrix_main(const char* name, int argc, char** argv)
{
    Parameters params;
    cov_matrix_params(name, argc, argv, params);

    string ifn_libraries = params.get_string("ifn_libraries");
    string ifn_segments = params.get_string("ifn_segments");
    string ifn_fasta = params.get_string("ifn_fasta");
    string ofn_mat = params.get_string("ofn_mat");
    string ofn_fasta = params.get_string("ofn_fasta");
    bool actual_nts = params.get_bool("actual_nts");
    bool should_create_fasta = params.get_bool("should_create_fasta");
    int min_segment_length = params.get_int("min_segment_length");
    
    ClipMode clip_mode = string_to_clip_mode(params.get_string("clip_mode"));
    int clip_margin = params.get_int("clip_margin");
    double min_mutations_percent = params.get_double("min_mutations_percent");
    double max_mutations_percent = params.get_double("max_mutations_percent");
    int min_alignment_length = params.get_int("min_alignment_length");
    int max_alignment_length = params.get_int("max_alignment_length");
    int min_indel_length = params.get_int("min_indel_length");
    string ofn_lib_map = params.get_string("ofn_lib_map");

    CovMatrix cov_matrix;
    cov_matrix.compute(ifn_libraries, ifn_segments, ifn_fasta, ofn_mat, ofn_fasta,
                      actual_nts, should_create_fasta,
                      min_segment_length, clip_mode, clip_margin,
                      min_mutations_percent, max_mutations_percent,
                      min_alignment_length, max_alignment_length, min_indel_length,
                      ofn_lib_map);

    return 0;
}

