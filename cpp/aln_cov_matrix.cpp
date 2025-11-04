#include "CovMatrix.h"
#include "Params.h"
#include <iostream>

using namespace std;

void cov_matrix_params(const char* name, int argc, char** argv, Parameters& params)
{
    params.add_parser("ifn_libraries", new ParserFilename("library table with lib_id and aln_fn columns"), true);
    params.add_parser("ifn_segments", new ParserFilename("segment table with segment, contig, start, end columns"), true);
    params.add_parser("ifn_fasta", new ParserFilename("contig fasta file"), true);
    params.add_parser("ofn_mat", new ParserFilename("output coverage matrix file"), true);
    params.add_parser("ofn_fasta", new ParserFilename("output segment fasta file"), true);

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

    CovMatrix cov_matrix;
    cov_matrix.compute(ifn_libraries, ifn_segments, ifn_fasta, ofn_mat, ofn_fasta);

    return 0;
}

