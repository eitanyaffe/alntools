#include "Params.h"
#include "SegmentationManager.h"
#include "alignment_store.h"
#include "utils.h"
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

void segments_params(const char* name, int argc, char** argv, Parameters& params)
{
    params.add_parser("ifn_aln", new ParserFilename("input ALN file (single library mode)"), false);
    params.add_parser("ifn_libraries", new ParserFilename("input table with library definitions (multi-library mode)"), false);
    params.add_parser("ifn_contig_table", new ParserFilename("contig table with contig and length columns"), true);
    params.add_parser("odir", new ParserFilename("output directory"), false);
    params.add_parser("max_margin", new ParserInteger("maximum margin tolerance for breakpoint clustering (default 20)", 20), false);
    params.add_parser("min_anchor_length", new ParserInteger("minimum anchor alignment length (default 1000)", 1000), false);
    params.add_parser("min_dangle_length", new ParserInteger("minimum dangle alignment length (default 1000)", 1000), false);
    params.add_parser("max_anchor_mutations_percent", new ParserDouble("maximum mutations percentage for anchor alignments (default 0.1)", 0.1), false);
    params.add_parser("min_alignment_distance", new ParserInteger("minimum distance between anchor and dangle on same contig (default 200)", 200), false);
    params.add_parser("min_breakpoint_read_support", new ParserInteger("minimum read support for selecting breakpoints (default 2)", 2), false);
    params.add_parser("min_breakpoint_frequency", new ParserDouble("minimum frequency for selecting breakpoints (default 0.2)", 0.2), false);
    params.add_parser("min_segment_length", new ParserInteger("minimum segment length for filtering breakpoints (default 200)", 200), false);
    params.add_parser("max_segment_length", new ParserInteger("maximum segment length before splitting, 0 to disable (default 500000)", 500000), false);

    if (argc == 1) {
        params.usage(name);
        exit(1);
    }

    // read command line params
    params.read(argc, argv);
    params.parse();
    
    // validate that exactly one of ifn_aln or ifn_libraries is provided
    string ifn_aln = params.get_string("ifn_aln");
    string ifn_libraries = params.get_string("ifn_libraries");
    
    if (ifn_aln.empty() && ifn_libraries.empty()) {
        cerr << "error: either ifn_aln (single library) or ifn_libraries (multi-library) must be provided" << endl;
        exit(1);
    }
    
    if (!ifn_aln.empty() && !ifn_libraries.empty()) {
        cerr << "error: provide either ifn_aln OR ifn_libraries, not both" << endl;
        exit(1);
    }
    
    params.verify_mandatory();

    // validate parameters
    int max_margin = params.get_int("max_margin");
    int min_anchor_length = params.get_int("min_anchor_length");
    int min_dangle_length = params.get_int("min_dangle_length");
    double max_anchor_mutations_percent = params.get_double("max_anchor_mutations_percent");
    int min_alignment_distance = params.get_int("min_alignment_distance");
    int min_breakpoint_read_support = params.get_int("min_breakpoint_read_support");
    double min_breakpoint_frequency = params.get_double("min_breakpoint_frequency");
    int min_segment_length = params.get_int("min_segment_length");
    int max_segment_length = params.get_int("max_segment_length");

    if (max_margin < 0) {
        cerr << "error: max_margin must be non-negative" << endl;
        exit(1);
    }

    if (min_anchor_length <= 0) {
        cerr << "error: min_anchor_length must be positive" << endl;
        exit(1);
    }

    if (min_dangle_length <= 0) {
        cerr << "error: min_dangle_length must be positive" << endl;
        exit(1);
    }

    if (max_anchor_mutations_percent < 0.0 || max_anchor_mutations_percent > 100.0) {
        cerr << "error: max_anchor_mutations_percent must be between 0.0 and 100.0" << endl;
        exit(1);
    }

    if (min_alignment_distance < 0) {
        cerr << "error: min_alignment_distance must be non-negative" << endl;
        exit(1);
    }

    if (min_breakpoint_read_support < 0) {
        cerr << "error: min_breakpoint_read_support must be non-negative" << endl;
        exit(1);
    }

    if (min_breakpoint_frequency < 0.0 || min_breakpoint_frequency > 1.0) {
        cerr << "error: min_breakpoint_frequency must be between 0.0 and 1.0" << endl;
        exit(1);
    }

    if (min_segment_length <= 0) {
        cerr << "error: min_segment_length must be positive" << endl;
        exit(1);
    }

    if (max_segment_length < 0) {
        cerr << "error: max_segment_length must be non-negative (0 to disable splitting)" << endl;
        exit(1);
    }

    if (max_segment_length > 0 && max_segment_length < min_segment_length) {
        cerr << "error: max_segment_length must be >= min_segment_length or 0 to disable" << endl;
        exit(1);
    }

    params.print(cout);
}

int segments_main(const char* name, int argc, char** argv)
{
    Parameters params;
    segments_params(name, argc, argv, params);

    string ifn_aln = params.get_string("ifn_aln");
    string ifn_libraries = params.get_string("ifn_libraries");
    string ifn_contig_table = params.get_string("ifn_contig_table");
    string odir = params.get_string("odir");
    int max_margin = params.get_int("max_margin");
    int min_anchor_length = params.get_int("min_anchor_length");
    int min_dangle_length = params.get_int("min_dangle_length");
    double max_anchor_mutations_percent = params.get_double("max_anchor_mutations_percent");
    int min_alignment_distance = params.get_int("min_alignment_distance");
    int min_breakpoint_read_support = params.get_int("min_breakpoint_read_support");
    double min_breakpoint_frequency = params.get_double("min_breakpoint_frequency");
    int min_segment_length = params.get_int("min_segment_length");
    int max_segment_length = params.get_int("max_segment_length");

    cout << "segments command called:" << endl;
    if (!ifn_aln.empty()) {
        cout << "  mode: single library" << endl;
        cout << "  ifn_aln: " << ifn_aln << endl;
    } else {
        cout << "  mode: multi-library" << endl;
        cout << "  ifn_libraries: " << ifn_libraries << endl;
    }
    cout << "  ifn_contig_table: " << ifn_contig_table << endl;
    cout << "  odir: " << odir << endl;
    cout << "  max_margin: " << max_margin << endl;
    cout << "  min_anchor_length: " << min_anchor_length << endl;
    cout << "  min_dangle_length: " << min_dangle_length << endl;
    cout << "  max_anchor_mutations_percent: " << max_anchor_mutations_percent << endl;
    cout << "  min_alignment_distance: " << min_alignment_distance << endl;
    cout << "  min_breakpoint_read_support: " << min_breakpoint_read_support << endl;
    cout << "  min_breakpoint_frequency: " << min_breakpoint_frequency << endl;
    cout << "  min_segment_length: " << min_segment_length << endl;
    cout << "  max_segment_length: " << max_segment_length << endl;

    // create library files map
    map<string, string> library_files;
    
    if (!ifn_aln.empty()) {
        // single file mode
        library_files["sample"] = ifn_aln;
        cout << "single library mode: using " << ifn_aln << " as library 'sample'" << endl;
    } else {
        // multi-library mode
        cout << "loading libraries from " << ifn_libraries << endl;
        vector<string> library_ids;  // not used but required by function signature
        read_library_table(ifn_libraries, library_files, library_ids);
    }
    
    // load all ALN stores
    map<string, AlignmentStore> stores;
    for (const auto& entry : library_files) {
        const string& lib_id = entry.first;
        const string& aln_file = entry.second;
        
        cout << "loading library " << lib_id << " from " << aln_file << endl;
        AlignmentStore store;
        store.load(aln_file);
        
        // build read-to-alignment index
        if (!store.is_read_alignment_index_built()) {
            store.init_read_alignment_index();
        }
        
        stores[lib_id] = std::move(store);
    }
    
    // create segmentation manager
    SegmentationManager manager(stores, ifn_contig_table, max_margin, min_anchor_length,
                                min_dangle_length, max_anchor_mutations_percent,
                                min_alignment_distance, min_breakpoint_read_support,
                                min_breakpoint_frequency, min_segment_length, max_segment_length);
    
    // execute segmentation analysis
    if (odir.empty()) {
        cerr << "error: odir must be provided" << endl;
        exit(1);
    }
    manager.execute();
    manager.write_to_csv(odir);

    cout << "segments command completed successfully" << endl;
    return 0;
}

