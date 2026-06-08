#include "GetLocalDeletions.h"
#include "Params.h"
#include <iostream>

using namespace std;

int get_local_deletions_main(const char* name, int argc, char** argv)
{
    Parameters params;
    params.add_parser("ifn_aln",          new ParserFilename("input ALN file"), true);
    params.add_parser("ifn_intervals",    new ParserFilename("input intervals file (contig, start, end)"), true);
    params.add_parser("ofn",             new ParserFilename("output tab-delimited file"), true);
    params.add_parser("margin",          new ParserInteger("max coordinate difference allowed between deletion and query interval endpoints (default 10)", 10), false);
    params.add_parser("flank",           new ParserInteger("window in bp around deletion for reporting nearby mutations (default 100)", 100), false);
    params.add_parser("min_indel_length", new ParserInteger("minimum indel length to include in flanking mutation reporting (default 3)", 3), false);

    if (argc == 1) {
        params.usage(name);
        exit(1);
    }

    params.read(argc, argv);
    params.parse();
    params.verify_mandatory();
    params.print(cout);

    string ifn_aln          = params.get_string("ifn_aln");
    string ifn_intervals    = params.get_string("ifn_intervals");
    string ofn              = params.get_string("ofn");
    int    margin           = params.get_int("margin");
    int    flank            = params.get_int("flank");
    int    min_indel_length = params.get_int("min_indel_length");

    GetLocalDeletions gld;
    gld.compute(ifn_aln, ifn_intervals, ofn, margin, flank, min_indel_length);

    return 0;
}
