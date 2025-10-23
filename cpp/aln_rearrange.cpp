#include "Params.h"
#include "Rearrange.h"
#include "RearrangeVerify.h"
#include "Sequences.h"
#include "alignment_store.h"
#include "utils.h"
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace std;

void rearrange_params(const char* name, int argc, char** argv, Parameters& params)
{
    params.add_parser("ifn_aln", new ParserFilename("input ALN file (single library mode)"), false);
    params.add_parser("ifn_libraries", new ParserFilename("input table with library definitions (multi-library mode)"), false);
    params.add_parser("ifn_intervals", new ParserFilename("input table with query contig intervals"), false);
    params.add_parser("ofn_prefix", new ParserFilename("output prefix"), false);
    params.add_parser("max_margin", new ParserInteger("maximum margin tolerance for all geometric constraints (default 10)", 10), false);
    params.add_parser("min_element_length", new ParserInteger("minimum element length for deletions (default 50)", 50), false);
    params.add_parser("min_anchor_length", new ParserInteger("minimum anchor alignment length (default 200)", 200), false);
    params.add_parser("max_anchor_mutations_percent", new ParserDouble("maximum mutations percentage for anchor alignments (default 1.0)", 1.0), false);
    params.add_parser("max_element_mutation_percent", new ParserDouble("maximum mutations percentage for element alignments (default 1.0)", 1.0), false);
    params.add_parser("min_indel_length", new ParserInteger("minimum indel length to include in mutation density calculations (default 3)", 3), false);
    params.add_parser("should_verify", new ParserBoolean("verify rearrangement events against reference sequences (default false)", false), false);
    params.add_parser("verify_error_mode", new ParserString("verification error handling: exit_on_error, exit_if_error, warning_only (default exit_on_error)", "exit_on_error"), false);
    params.add_parser("ifn_contigs", new ParserFilename("input contig FASTA file (required if should_verify=true or resolve_seams requires reference)"), false);
    params.add_parser("resolve_seams", new ParserString("seam resolution mode: no, reads_only, reference_only, complete (default no)", "no"), false);

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
    int min_element_length = params.get_int("min_element_length");
    int min_anchor_length = params.get_int("min_anchor_length");
    double max_anchor_mutations_percent = params.get_double("max_anchor_mutations_percent");
    double max_element_mutation_percent = params.get_double("max_element_mutation_percent");
    int min_indel_length = params.get_int("min_indel_length");

    if (max_margin < 0) {
        cerr << "error: max_margin must be non-negative" << endl;
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

    if (max_anchor_mutations_percent < 0.0 || max_anchor_mutations_percent > 100.0) {
        cerr << "error: max_anchor_mutations_percent must be between 0.0 and 100.0" << endl;
        exit(1);
    }

    if (max_element_mutation_percent < 0.0 || max_element_mutation_percent > 100.0) {
        cerr << "error: max_element_mutation_percent must be between 0.0 and 100.0" << endl;
        exit(1);
    }

    if (min_indel_length <= 0) {
        cerr << "error: min_indel_length must be positive" << endl;
        exit(1);
    }

    bool should_verify = params.get_bool("should_verify");
    string ifn_contigs = params.get_string("ifn_contigs");
    
    // basic validation - detailed validation will be done in main function
    if (should_verify && ifn_contigs.empty()) {
        cerr << "error: ifn_contigs is required when should_verify=true" << endl;
        exit(1);
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
    int max_margin = params.get_int("max_margin");
    int min_element_length = params.get_int("min_element_length");
    int min_anchor_length = params.get_int("min_anchor_length");
    double max_anchor_mutations_percent = params.get_double("max_anchor_mutations_percent");
    double max_element_mutation_percent = params.get_double("max_element_mutation_percent");
    int min_indel_length = params.get_int("min_indel_length");
    bool should_verify = params.get_bool("should_verify");
    string resolve_seams_str = params.get_string("resolve_seams");
    ResolveSeams resolve_seams = string_to_resolve_seams(resolve_seams_str);
    string ifn_contigs = params.get_string("ifn_contigs");
    string verify_error_mode_str = params.get_string("verify_error_mode");
    VerifyErrorMode verify_error_mode = string_to_verify_error_mode(verify_error_mode_str);
    
    // validation based on resolve_seams mode
    if (should_verify || resolve_seams == ResolveSeams::REFERENCE_ONLY || resolve_seams == ResolveSeams::COMPLETE) {
        if (ifn_contigs.empty()) {
            cerr << "error: ifn_contigs is required for verification or reference seam resolution" << endl;
            exit(1);
        }
    }

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
    cout << "  max_margin: " << max_margin << endl;
    cout << "  min_element_length: " << min_element_length << endl;
    cout << "  min_anchor_length: " << min_anchor_length << endl;
    cout << "  max_anchor_mutations_percent: " << max_anchor_mutations_percent << endl;
    cout << "  max_element_mutation_percent: " << max_element_mutation_percent << endl;
    cout << "  min_indel_length: " << min_indel_length << endl;
    cout << "  should_verify: " << (should_verify ? "true" : "false") << endl;
    if (should_verify) {
        cout << "  verify_error_mode: " << verify_error_mode_to_string(verify_error_mode) << endl;
    }
    cout << "  resolve_seams: " << resolve_seams_to_string(resolve_seams) << endl;
    if (should_verify || resolve_seams == ResolveSeams::REFERENCE_ONLY || resolve_seams == ResolveSeams::COMPLETE) {
        cout << "  ifn_contigs: " << ifn_contigs << endl;
    }

    // Load intervals (optional)
    vector<Interval> intervals;
    if (!ifn_intervals.empty()) {
        read_intervals(ifn_intervals, intervals);
        cout << "found " << intervals.size() << " intervals in " << ifn_intervals << endl;
    } else {
        cout << "no intervals file provided - processing all reads" << endl;
    }

    // Always use Rearrange - create single library for ifn_aln mode
    map<string, string> library_files;
    map<string, string> library_read_files;  // lib_id -> read_file
    
    if (!ifn_aln.empty()) {
        // Single file mode - create single library with id "sample"
        library_files["sample"] = ifn_aln;
        // Note: read files will be handled by library table for multi-library mode
        // For single library mode, read sequences would need to be loaded separately
        cout << "single library mode: using " << ifn_aln << " as library 'sample'" << endl;
    } else {
        // Multi-library mode - read libraries table with extended format
        cout << "loading libraries from " << ifn_libraries << endl;
        map<string, LibraryInfo> library_info = read_library_table_extended(ifn_libraries);
        
        for (const auto& entry : library_info) {
            const string& lib_id = entry.first;
            const LibraryInfo& info = entry.second;
            
            library_files[lib_id] = info.aln_file;
            if (!info.read_file.empty()) {
                library_read_files[lib_id] = info.read_file;
            }
        }
    }
    
    // load all ALN stores
    map<string, AlignmentStore> stores;
    for (const auto& entry : library_files) {
        const string& lib_id = entry.first;
        const string& aln_file = entry.second;
        
        cout << "loading library " << lib_id << " from " << aln_file << endl;
        AlignmentStore store;
        store.load(aln_file);
        
        // count short indels for this store
        store.count_short_indels(min_indel_length);
        
        // build read-to-alignment index
        if (!store.is_read_alignment_index_built()) {
            store.init_read_alignment_index();
        }
        
        stores[lib_id] = std::move(store);
    }
    
    // create sequence objects and verifier if needed
    unique_ptr<AssemblySequences> assembly_sequences;
    unique_ptr<ReadSequences> read_sequences;
    RearrangeVerify* verifier = nullptr;
    unique_ptr<RearrangeVerify> verify_obj;
    
    if (should_verify || resolve_seams == ResolveSeams::REFERENCE_ONLY || resolve_seams == ResolveSeams::COMPLETE) {
        cout << "loading assembly sequences..." << endl;
        assembly_sequences = make_unique<AssemblySequences>();
        assembly_sequences->load_from_file(ifn_contigs);
    }
    
    if (should_verify || resolve_seams == ResolveSeams::READS_ONLY || resolve_seams == ResolveSeams::COMPLETE) {
        cout << "registering read sequence files..." << endl;
        read_sequences = make_unique<ReadSequences>();
        // register all library files (no loading into memory)
        for (const auto& entry : library_read_files) {
            const string& lib_id = entry.first;
            const string& read_file = entry.second;
            read_sequences->register_lib_file(lib_id, read_file);
        }
    }
    
    if (should_verify) {
        cout << "verification enabled" << endl;
        verify_obj = make_unique<RearrangeVerify>(assembly_sequences.get(), verify_error_mode);
        verifier = verify_obj.get();
    }
    
    // create rearrangement manager
    Rearrange manager(stores, intervals, verifier, resolve_seams, assembly_sequences.get(), read_sequences.get(), 
                     max_margin, min_element_length, min_anchor_length, max_anchor_mutations_percent, max_element_mutation_percent, true, ofn_prefix);
    
    // Note: read sequences are now loaded directly into sequence objects above
    
    // execute rearrangement analysis
    manager.execute();
    manager.write_to_csv(ofn_prefix);

    cout << "rearrange command completed successfully" << endl;
    return 0;
}
