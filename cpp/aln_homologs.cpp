#include "Homologs.h"
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

void homologs_usage(const char* name) {
    fprintf(stderr, "usage: %s [options]\n", name);
    fprintf(stderr, "options:\n");
    fprintf(stderr, "  -ifn_fasta <filename>    input assembly FASTA file\n");
    fprintf(stderr, "  -contig <string>         query contig name\n");
    fprintf(stderr, "  -start <int>             query start position (1-based)\n");
    fprintf(stderr, "  -end <int>               query end position (1-based, inclusive)\n");
    fprintf(stderr, "  -k <int>                 kmer size (default: 21)\n");
    fprintf(stderr, "  -num_kmers <int>         number of kmers to extract (default: 10)\n");
    fprintf(stderr, "  -threshold <double>      minimum percentage of kmers required (default: 80.0)\n");
    fprintf(stderr, "  -num_threads <int>       number of threads (0 for auto, default: 0)\n");
    fprintf(stderr, "  -ofn <filename>          output regions file\n");
}

int homologs_main(const char* name, int argc, char** argv) {
    // default parameters
    string ifn_fasta = "";
    string query_contig = "";
    uint32_t query_start = 0;
    uint32_t query_end = 0;
    uint32_t k = 21;
    uint32_t num_kmers = 10;
    double threshold = 80.0;
    int num_threads = 0;
    string ofn = "";
    
    // parse command line arguments
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "-ifn_fasta" && i + 1 < argc) {
            ifn_fasta = argv[++i];
        } else if (arg == "-contig" && i + 1 < argc) {
            query_contig = argv[++i];
        } else if (arg == "-start" && i + 1 < argc) {
            query_start = atoi(argv[++i]);
        } else if (arg == "-end" && i + 1 < argc) {
            query_end = atoi(argv[++i]);
        } else if (arg == "-k" && i + 1 < argc) {
            k = atoi(argv[++i]);
        } else if (arg == "-num_kmers" && i + 1 < argc) {
            num_kmers = atoi(argv[++i]);
        } else if (arg == "-threshold" && i + 1 < argc) {
            threshold = atof(argv[++i]);
        } else if (arg == "-num_threads" && i + 1 < argc) {
            num_threads = atoi(argv[++i]);
        } else if (arg == "-ofn" && i + 1 < argc) {
            ofn = argv[++i];
        } else {
            fprintf(stderr, "error: unknown argument: %s\n", arg.c_str());
            homologs_usage(name);
            return 1;
        }
    }
    
    // validate required parameters
    if (ifn_fasta.empty()) {
        fprintf(stderr, "error: -ifn_fasta is required\n");
        homologs_usage(name);
        return 1;
    }
    
    if (query_contig.empty()) {
        fprintf(stderr, "error: -contig is required\n");
        homologs_usage(name);
        return 1;
    }
    
    if (query_start == 0 || query_end == 0) {
        fprintf(stderr, "error: -start and -end are required\n");
        homologs_usage(name);
        return 1;
    }
    
    if (query_start > query_end) {
        fprintf(stderr, "error: start position must be <= end position\n");
        return 1;
    }
    
    if (ofn.empty()) {
        fprintf(stderr, "error: -ofn is required\n");
        homologs_usage(name);
        return 1;
    }
    
    if (k == 0 || num_kmers == 0) {
        fprintf(stderr, "error: k and num_kmers must be > 0\n");
        return 1;
    }
    
    if (threshold <= 0.0 || threshold > 100.0) {
        fprintf(stderr, "error: threshold must be between 0 and 100\n");
        return 1;
    }
    
    try {
        cout << "loading assembly sequences from " << ifn_fasta << endl;
        cout << "extracting " << num_kmers << " kmers of size " << k 
             << " from " << query_contig << ":" << query_start << "-" << query_end << endl;
        cout << "searching for regions with >= " << threshold << "% kmer coverage" << endl;
        
        // create Homologs object with threading support
        Homologs homologs(num_threads);
        
        vector<ContigRegion> regions = homologs.search_homologs(ifn_fasta, query_contig, query_start, query_end, 
                                                               k, num_kmers, threshold);
        
        cout << "found " << regions.size() << " homologous regions" << endl;
        
        homologs.write_results(regions, ofn);
        cout << "results written to " << ofn << endl;
        
    } catch (const exception& e) {
        fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    
    return 0;
}
