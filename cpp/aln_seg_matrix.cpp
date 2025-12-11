#include "SegMatrix.h"
#include "alignment_store.h"
#include "Params.h"
#include "utils.h"
#include <iostream>
#include <filesystem>

using namespace std;

void seg_matrix_params(const char* name, int argc, char** argv, Parameters& params)
{
    params.add_parser("ifn_aln", new ParserFilename("input ALN file"), true);
    params.add_parser("ifn_segments", new ParserFilename("segment table with segment, contig, start, end columns"), true);
    params.add_parser("odir", new ParserFilename("output directory"), true);
    params.add_parser("max_mutation_percent", new ParserDouble("maximum mutation percentage for side coverage (default 1.0)", 1.0), false);
    params.add_parser("max_adjacency_distance", new ParserInteger("maximum distance for adjacency in bp (default 200)", 200), false);
    params.add_parser("max_margin", new ParserInteger("maximum overlap between alignments to skip in bp (default 10)", 10), false);
    params.add_parser("min_indel_length", new ParserInteger("minimum indel length to include in mutation density calculations (default 3)", 3), false);
    params.add_parser("side_length", new ParserInteger("side length in bp for segment coverage analysis (default 1000)", 1000), false);
    params.add_parser("side_margin", new ParserInteger("margin in bp between segment edge and side region (default 200)", 200), false);
    params.add_parser("output_read_details", new ParserBoolean("output detailed read association files (default F)", false), false);
    params.add_parser("ifn_segment_clusters", new ParserFilename("optional segment-cluster mapping table with segment, segment_cluster, strand columns"), false);

    if (argc == 1) {
        params.usage(name);
        exit(1);
    }

    params.read(argc, argv);
    params.parse();
    params.verify_mandatory();
    params.print(cout);
}

int seg_matrix_main(const char* name, int argc, char** argv)
{
    Parameters params;
    seg_matrix_params(name, argc, argv, params);

    string ifn_aln = params.get_string("ifn_aln");
    string ifn_segments = params.get_string("ifn_segments");
    string odir = params.get_string("odir");
    double max_mutation_percent = params.get_double("max_mutation_percent");
    int max_adjacency_distance = params.get_int("max_adjacency_distance");
    int max_margin = params.get_int("max_margin");
    int min_indel_length = params.get_int("min_indel_length");
    int side_length = params.get_int("side_length");
    int side_margin = params.get_int("side_margin");
    bool output_read_details = params.get_bool("output_read_details");
    string ifn_segment_clusters = params.is_defined("ifn_segment_clusters") ? params.get_string("ifn_segment_clusters") : "";

    // create output directory
    {
        std::error_code ec;
        std::filesystem::create_directories(odir, ec);
        if (ec) {
            cerr << "error: could not create output directory " << odir << ": " << ec.message() << endl;
            exit(1);
        }
    }

    // load alignment store
    cout << "loading alignment store: " << ifn_aln << endl;
    AlignmentStore store;
    store.load(ifn_aln);
    cout << "loaded " << store.get_read_count() << " reads and " 
         << store.get_alignment_count() << " alignments" << endl;

    SegMatrix seg_matrix(store);
    seg_matrix.compute(ifn_segments, odir,
                      max_mutation_percent, max_adjacency_distance,
                      max_margin, min_indel_length,
                      side_length, side_margin, output_read_details,
                      ifn_segment_clusters);

    cout << "seg_matrix command completed successfully" << endl;
    return 0;
}

