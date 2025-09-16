#include "Params.h"
#include "Rearrange.h"
#include "RearrangeManager.h"
#include "VerifyRearrange.h"
#include "alignment_store.h"
#include "utils.h"
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

void rearrange_params(const char* name, int argc, char** argv, Parameters& params)
{
    params.add_parser("ifn_aln", new ParserFilename("input ALN file (single library mode)"), false);
    params.add_parser("ifn_libraries", new ParserFilename("input table with library definitions (multi-library mode)"), false);
    params.add_parser("ifn_intervals", new ParserFilename("input table with query contig intervals"), false);
    params.add_parser("ofn_prefix", new ParserFilename("output prefix"), false);
    params.add_parser("max_gap", new ParserInteger("maximum gap tolerance for all geometric constraints (default 10)", 10), false);
    params.add_parser("min_element_length", new ParserInteger("minimum element length for deletions (default 50)", 50), false);
    params.add_parser("min_anchor_length", new ParserInteger("minimum anchor alignment length (default 200)", 200), false);
    params.add_parser("max_mutations_percent", new ParserDouble("maximum mutations percentage for all alignments (default 0.01)", 0.01), false);
    params.add_parser("should_verify", new ParserBoolean("verify rearrangement events against reference sequences (default false)", false), false);
    params.add_parser("ifn_contigs", new ParserFilename("input contig FASTA file (required if should_verify=true)"), false);
    params.add_parser("ifn_reads", new ParserFilename("input read FASTQ file (required if should_verify=true)"), false);

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
    int max_gap = params.get_int("max_gap");
    int min_element_length = params.get_int("min_element_length");
    int min_anchor_length = params.get_int("min_anchor_length");
    double max_mutations_percent = params.get_double("max_mutations_percent");

    if (max_gap < 0) {
        cerr << "error: max_gap must be non-negative" << endl;
        exit(1);
    }

    if (min_element_length <= 0) {
        cerr << "error: min_element_length must be positive" << endl;
        exit(1);
    }

    if (min_anchor_length <= 0) {
        cerr << "error: min_anchor_length must be positive" << endl;
        exit(1);
    }

    if (max_mutations_percent < 0.0 || max_mutations_percent > 1.0) {
        cerr << "error: max_mutations_percent must be between 0.0 and 1.0" << endl;
        exit(1);
    }

    bool should_verify = params.get_bool("should_verify");
    string ifn_contigs = params.get_string("ifn_contigs");
    string ifn_reads = params.get_string("ifn_reads");
    
    if (should_verify) {
        if (ifn_contigs.empty()) {
            cerr << "error: ifn_contigs is required when should_verify=true" << endl;
            exit(1);
        }
        if (ifn_reads.empty()) {
            cerr << "error: ifn_reads is required when should_verify=true" << endl;
            exit(1);
        }
    }

    params.print(cout);
}

int rearrange_main(const char* name, int argc, char** argv)
{
    Parameters params;
    rearrange_params(name, argc, argv, params);

    string ifn_aln = params.get_string("ifn_aln");
    string ifn_libraries = params.get_string("ifn_libraries");
    string ifn_intervals = params.get_string("ifn_intervals");
    string ofn_prefix = params.get_string("ofn_prefix");
    int max_gap = params.get_int("max_gap");
    int min_element_length = params.get_int("min_element_length");
    int min_anchor_length = params.get_int("min_anchor_length");
    double max_mutations_percent = params.get_double("max_mutations_percent");
    bool should_verify = params.get_bool("should_verify");
    string ifn_contigs = params.get_string("ifn_contigs");
    string ifn_reads = params.get_string("ifn_reads");

    cout << "rearrange command called:" << endl;
    if (!ifn_aln.empty()) {
        cout << "  mode: single library" << endl;
        cout << "  ifn_aln: " << ifn_aln << endl;
    } else {
        cout << "  mode: multi-library" << endl;
        cout << "  ifn_libraries: " << ifn_libraries << endl;
    }
    cout << "  ifn_intervals: " << ifn_intervals << endl;
    cout << "  ofn_prefix: " << ofn_prefix << endl;
    cout << "  max_gap: " << max_gap << endl;
    cout << "  min_element_length: " << min_element_length << endl;
    cout << "  min_anchor_length: " << min_anchor_length << endl;
    cout << "  max_mutations_percent: " << max_mutations_percent << endl;
    cout << "  should_verify: " << (should_verify ? "true" : "false") << endl;
    if (should_verify) {
        cout << "  ifn_contigs: " << ifn_contigs << endl;
        cout << "  ifn_reads: " << ifn_reads << endl;
    }

    // Load intervals (optional)
    vector<Interval> intervals;
    if (!ifn_intervals.empty()) {
        read_intervals(ifn_intervals, intervals);
        cout << "found " << intervals.size() << " intervals in " << ifn_intervals << endl;
    } else {
        cout << "no intervals file provided - processing all reads" << endl;
    }

    // Always use RearrangeManager - create single library for ifn_aln mode
    map<string, string> library_files;
    
    if (!ifn_aln.empty()) {
        // Single file mode - create single library with id "sample"
        library_files["sample"] = ifn_aln;
        cout << "single library mode: using " << ifn_aln << " as library 'sample'" << endl;
    } else {
        // Multi-library mode - read libraries table
        cout << "loading libraries from " << ifn_libraries << endl;
        library_files = read_library_table(ifn_libraries);
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
            cout << "building read-to-alignment index for library " << lib_id << endl;
            store.init_read_alignment_index();
        }
        
        stores[lib_id] = std::move(store);
    }
    
    // create verifier if needed
    VerifyRearrange* verifier = nullptr;
    VerifyRearrange verify_obj;
    if (should_verify) {
        cout << "verification enabled - loading sequences..." << endl;
        verify_obj.load_sequences(ifn_contigs, ifn_reads);
        verifier = &verify_obj;
    }
    
    // create and execute rearrangement manager
    RearrangeManager manager(stores, intervals, verifier, max_gap, min_element_length, min_anchor_length, max_mutations_percent);
    manager.execute();
    manager.write_to_csv(ofn_prefix);

    cout << "rearrange command completed successfully" << endl;
    return 0;
}
