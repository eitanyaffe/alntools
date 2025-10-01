#include "Homologs.h"
#include "Sequences.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

Homologs::Homologs() {}

Homologs::~Homologs() {}

vector<string> Homologs::extract_query_kmers(const string& sequence, uint32_t start, uint32_t end, 
                                            uint32_t k, uint32_t num_kmers) {
    vector<string> kmers;
    
    if (end <= start || end - start + 1 < k) {
        return kmers;
    }
    
    uint32_t region_length = end - start + 1;
    uint32_t max_possible_kmers = region_length - k + 1;
    
    if (max_possible_kmers == 0) {
        return kmers;
    }
    
    // adjust num_kmers if region is too small
    uint32_t actual_num_kmers = min(num_kmers, max_possible_kmers);
    
    if (actual_num_kmers == 1) {
        // single kmer at the beginning
        string kmer = sequence.substr(start - 1, k);
        kmers.push_back(kmer);
    } else {
        // evenly space kmers across the region
        double step = (double)(max_possible_kmers - 1) / (actual_num_kmers - 1);
        
        for (uint32_t i = 0; i < actual_num_kmers; i++) {
            uint32_t pos = start - 1 + (uint32_t)(i * step);
            string kmer = sequence.substr(pos, k);
            kmers.push_back(kmer);
        }
    }
    
    return kmers;
}

unordered_map<string, vector<KmerHit>> Homologs::build_kmer_index(const AssemblySequences& assembly, uint32_t k) {
    unordered_map<string, vector<KmerHit>> kmer_index;
    
    auto contig_ids = assembly.get_contig_ids();
    
    for (const string& contig_id : contig_ids) {
        string sequence = assembly.get_sequence(contig_id);
        
        if (sequence.length() < k) {
            continue;
        }
        
        // extract all kmers from this contig
        for (uint32_t pos = 0; pos <= sequence.length() - k; pos++) {
            string kmer = sequence.substr(pos, k);
            
            // convert to uppercase for consistency
            transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
            
            kmer_index[kmer].push_back({contig_id, pos + 1}); // 1-based position
        }
    }
    
    return kmer_index;
}

vector<ContigRegion> Homologs::find_homologous_regions(const vector<string>& query_kmers,
                                                      const unordered_map<string, vector<KmerHit>>& kmer_index,
                                                      double threshold_percent) {
    
    // collect all hits per contig
    unordered_map<string, vector<uint32_t>> contig_hits;
    
    uint32_t total_query_kmers = query_kmers.size();
    uint32_t min_required_kmers = (uint32_t)ceil(total_query_kmers * threshold_percent / 100.0);
    
    for (const string& kmer : query_kmers) {
        string upper_kmer = kmer;
        transform(upper_kmer.begin(), upper_kmer.end(), upper_kmer.begin(), ::toupper);
        
        auto it = kmer_index.find(upper_kmer);
        if (it != kmer_index.end()) {
            for (const KmerHit& hit : it->second) {
                contig_hits[hit.contig].push_back(hit.position);
            }
        }
    }
    
    vector<ContigRegion> regions;
    
    // process each contig with hits
    for (auto& pair : contig_hits) {
        const string& contig = pair.first;
        vector<uint32_t>& positions = pair.second;
        
        if (positions.size() < min_required_kmers) {
            continue;
        }
        
        // sort positions
        sort(positions.begin(), positions.end());
        
        // find span from first to last hit
        uint32_t start = positions.front();
        uint32_t end = positions.back();
        uint32_t length = end - start + 1;
        uint32_t kmer_count = positions.size();
        double coverage_percent = (double)kmer_count / total_query_kmers * 100.0;
        
        regions.push_back({contig, start, end, length, kmer_count, coverage_percent});
    }
    
    // sort by coverage percentage (descending)
    sort(regions.begin(), regions.end(), 
         [](const ContigRegion& a, const ContigRegion& b) {
             return a.coverage_percent > b.coverage_percent;
         });
    
    return regions;
}

vector<ContigRegion> Homologs::search_homologs(const string& fasta_file,
                                              const string& query_contig,
                                              uint32_t query_start,
                                              uint32_t query_end,
                                              uint32_t k,
                                              uint32_t num_kmers,
                                              double threshold_percent) {
    
    // load assembly sequences
    AssemblySequences assembly;
    assembly.load_from_file(fasta_file);
    
    if (!assembly.has_contig(query_contig)) {
        throw runtime_error("contig '" + query_contig + "' not found in assembly");
    }
    
    // get query sequence
    string query_sequence = assembly.get_sequence(query_contig);
    if (query_end > query_sequence.length()) {
        throw runtime_error("end position " + to_string(query_end) + 
                          " exceeds contig length " + to_string(query_sequence.length()));
    }
    
    // extract query kmers
    vector<string> query_kmers = extract_query_kmers(query_sequence, query_start, query_end, k, num_kmers);
    
    if (query_kmers.empty()) {
        throw runtime_error("no kmers could be extracted from query region");
    }
    
    // build kmer index
    auto kmer_index = build_kmer_index(assembly, k);
    
    // find homologous regions
    return find_homologous_regions(query_kmers, kmer_index, threshold_percent);
}

void Homologs::write_results(const vector<ContigRegion>& regions, const string& output_file) {
    ofstream out(output_file);
    if (!out) {
        throw runtime_error("cannot open output file: " + output_file);
    }
    
    // write header
    out << "assembly\tcontig\tstart\tend\tdesc\tid" << endl;
    
    // write regions
    for (size_t i = 0; i < regions.size(); i++) {
        const ContigRegion& region = regions[i];
        
        // create description with length and coverage info
        string desc = "kmer_match_length_" + to_string(region.length) + 
                     "_coverage_" + to_string((int)region.coverage_percent) + "pct";
        
        string region_id = "homolog_" + to_string(i + 1);
        
        out << "assembly\t" << region.contig << "\t" << region.start << "\t" 
            << region.end << "\t" << desc << "\t" << region_id << endl;
    }
    
    out.close();
}

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
    fprintf(stderr, "  -ofn <filename>          output regions file\n");
}

