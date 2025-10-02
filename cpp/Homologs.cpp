#include "Homologs.h"
#include "Sequences.h"
#include "utils.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

Homologs::Homologs(int num_threads) : num_threads(num_threads) {
    // set number of threads for OpenMP
    if (num_threads <= 0) {
#ifdef _OPENMP
        this->num_threads = omp_get_max_threads();
#else
        this->num_threads = 1;
#endif
    }
    
#ifdef _OPENMP
    omp_set_num_threads(this->num_threads);
#endif
}

Homologs::~Homologs() {}

vector<string> Homologs::extract_query_kmers(const string& sequence, uint32_t start, uint32_t end, 
                                            uint32_t k, uint32_t num_kmers) {
    vector<string> kmers;
    
    cerr << "extracting query kmers from region " << start << "-" << end 
         << " (length: " << (end - start + 1) << ", k=" << k << ")" << endl;
    
    if (end <= start || end - start + 1 < k) {
        cerr << "region too small for kmer extraction" << endl;
        return kmers;
    }
    
    uint32_t region_length = end - start + 1;
    uint32_t max_possible_kmers = region_length - k + 1;
    
    if (max_possible_kmers == 0) {
        cerr << "no possible kmers in region" << endl;
        return kmers;
    }
    
    // adjust num_kmers if region is too small
    uint32_t actual_num_kmers = min(num_kmers, max_possible_kmers);
    cerr << "requested " << num_kmers << " kmers, extracting " << actual_num_kmers 
         << " (max possible: " << max_possible_kmers << ")" << endl;
    
    if (actual_num_kmers == 1) {
        // single kmer at the beginning
        string kmer = sequence.substr(start - 1, k);
        kmers.push_back(kmer);
        cerr << "extracted single kmer at position " << start << endl;
    } else {
        // evenly space kmers across the region
        double step = (double)(max_possible_kmers - 1) / (actual_num_kmers - 1);
        cerr << "spacing kmers with step size: " << step << endl;
        
        for (uint32_t i = 0; i < actual_num_kmers; i++) {
            uint32_t pos = start - 1 + (uint32_t)(i * step);
            string kmer = sequence.substr(pos, k);
            kmers.push_back(kmer);
        }
    }
    
    cerr << "successfully extracted " << kmers.size() << " query kmers" << endl;
    return kmers;
}

unordered_map<string, vector<uint32_t>> Homologs::build_contig_kmer_index(const string& sequence, uint32_t k) {
    unordered_map<string, vector<uint32_t>> kmer_index;
    
    if (sequence.length() < k) {
        return kmer_index;
    }
    
    // extract all kmers from this contig (both forward and reverse complement)
    for (uint32_t pos = 0; pos <= sequence.length() - k; pos++) {
        string kmer = sequence.substr(pos, k);
        
        // convert to uppercase for consistency
        transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
        
        // index forward kmer
        kmer_index[kmer].push_back(pos + 1); // 1-based position
        
        // index reverse complement kmer
        string rc_kmer = reverse_complement(kmer);
        if (rc_kmer != kmer) {  // avoid duplicates for palindromic kmers
            kmer_index[rc_kmer].push_back(pos + 1); // same position, different strand
        }
    }
    
    return kmer_index;
}

ContigRegion Homologs::process_contig(const string& contig_id,
                                     const string& sequence,
                                     const vector<string>& query_kmers,
                                     uint32_t k,
                                     double threshold_percent) {
    
    // build kmer index for this contig
    auto contig_kmer_index = build_contig_kmer_index(sequence, k);
    
    uint32_t total_query_kmers = query_kmers.size();
    uint32_t min_required_kmers = (uint32_t)ceil(total_query_kmers * threshold_percent / 100.0);
    
    // collect all hits in this contig
    vector<uint32_t> positions;
    uint32_t found_kmers = 0;
    
    for (const string& kmer : query_kmers) {
        string upper_kmer = kmer;
        transform(upper_kmer.begin(), upper_kmer.end(), upper_kmer.begin(), ::toupper);
        
        auto it = contig_kmer_index.find(upper_kmer);
        if (it != contig_kmer_index.end()) {
            found_kmers++;
            for (uint32_t pos : it->second) {
                positions.push_back(pos);
            }
        }
    }
    
    // return empty region if insufficient hits
    if (positions.size() < min_required_kmers) {
        return {contig_id, 0, 0, 0, 0, 0.0};
    }
    
    // sort positions and find span
    sort(positions.begin(), positions.end());
    
    uint32_t start = positions.front();
    uint32_t end = positions.back();
    uint32_t length = end - start + 1;
    uint32_t kmer_count = positions.size();
    double coverage_percent = (double)kmer_count / total_query_kmers * 100.0;
    
    return {contig_id, start, end, length, kmer_count, coverage_percent};
}

vector<ContigRegion> Homologs::search_homologs(const string& fasta_file,
                                              const string& query_contig,
                                              uint32_t query_start,
                                              uint32_t query_end,
                                              uint32_t k,
                                              uint32_t num_kmers,
                                              double threshold_percent) {
    
    cerr << "searching for homologs in contig " << query_contig 
         << " region " << query_start << "-" << query_end << " (from file: " << fasta_file << ")" << endl;
    
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
    
    // process contigs in parallel
    auto contig_ids = assembly.get_contig_ids();
    
    // filter out query contig and too-short contigs upfront
    vector<string> valid_contigs;
    for (const string& contig_id : contig_ids) {
        if (assembly.get_sequence(contig_id).length() >= k) {
            valid_contigs.push_back(contig_id);
        }
    }
    
#ifdef _OPENMP
    cerr << "processing " << valid_contigs.size() << " contigs using " << this->num_threads 
         << " threads (k=" << k << ", threshold=" << threshold_percent << "%)" << endl;
#else
    cerr << "processing " << valid_contigs.size() << " contigs using 1 thread (k=" << k 
         << ", threshold=" << threshold_percent << "%)" << endl;
#endif
    
    vector<ContigRegion> results;
    
#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        {
            cerr << "using " << omp_get_num_threads() << " threads for " << valid_contigs.size() << " contigs" << endl;
        }
        
        // thread-local results
        vector<ContigRegion> local_results;
        uint32_t local_processed = 0;
        
        #pragma omp for
        for (size_t i = 0; i < valid_contigs.size(); ++i) {
            const string& contig_id = valid_contigs[i];
            string sequence = assembly.get_sequence(contig_id);
            
            local_processed++;
            ContigRegion region = process_contig(contig_id, sequence, query_kmers, k, threshold_percent);
            
            // only keep regions that meet the threshold
            if (region.kmer_count > 0) {
                local_results.push_back(region);
            }
            
            // progress reporting (only from thread 0)
            if (omp_get_thread_num() == 0 && local_processed % 1000 == 0) {
                uint32_t estimated_total = local_processed * omp_get_num_threads();
                cerr << "progress: ~" << estimated_total << " contigs processed..." << endl;
            }
        }
        
        // merge local results into global results
        #pragma omp critical
        {
            results.insert(results.end(), local_results.begin(), local_results.end());
        }
    }
#else
    // fallback: sequential processing
    uint32_t processed_contigs = 0;
    
    for (const string& contig_id : valid_contigs) {
        string sequence = assembly.get_sequence(contig_id);
        
        processed_contigs++;
        ContigRegion region = process_contig(contig_id, sequence, query_kmers, k, threshold_percent);
        
        // only keep regions that meet the threshold
        if (region.kmer_count > 0) {
            results.push_back(region);
        }
        
        if (processed_contigs % 1000 == 0) {
            cerr << "processed " << processed_contigs << " contigs..." << endl;
        }
    }
#endif
    
    cerr << "processed " << valid_contigs.size() << " contigs, found " 
         << results.size() << " qualifying regions" << endl;
    
    // sort by coverage percentage (descending)
    sort(results.begin(), results.end(), 
         [](const ContigRegion& a, const ContigRegion& b) {
             return a.coverage_percent > b.coverage_percent;
         });
    
    if (!results.empty()) {
        cerr << "best match: contig " << results[0].contig 
             << " (" << results[0].coverage_percent << "% coverage, " 
             << results[0].length << " bp)" << endl;
    }
    
    cerr << "homolog search completed: " << results.size() << " regions found" << endl;
    return results;
}

vector<ContigRegion> Homologs::search_homologs(const AssemblySequences& assembly,
                                              const string& query_contig,
                                              uint32_t query_start,
                                              uint32_t query_end,
                                              uint32_t k,
                                              uint32_t num_kmers,
                                              double threshold_percent) {
    
    cerr << "searching for homologs in contig " << query_contig 
         << " region " << query_start << "-" << query_end << " (from loaded assembly)" << endl;
    
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
    
    // process contigs in parallel
    auto contig_ids = assembly.get_contig_ids();
    
    // filter out query contig and too-short contigs upfront
    vector<string> valid_contigs;
    for (const string& contig_id : contig_ids) {
        if (assembly.get_sequence(contig_id).length() >= k) {
            valid_contigs.push_back(contig_id);
        }
    }
    
#ifdef _OPENMP
    cerr << "processing " << valid_contigs.size() << " contigs using " << this->num_threads 
         << " threads (k=" << k << ", threshold=" << threshold_percent << "%)" << endl;
#else
    cerr << "processing " << valid_contigs.size() << " contigs using 1 thread (k=" << k 
         << ", threshold=" << threshold_percent << "%)" << endl;
#endif
    
    vector<ContigRegion> results;
    
#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        {
            cerr << "using " << omp_get_num_threads() << " threads for " << valid_contigs.size() << " contigs" << endl;
        }
        
        // thread-local results
        vector<ContigRegion> local_results;
        uint32_t local_processed = 0;
        
        #pragma omp for
        for (size_t i = 0; i < valid_contigs.size(); ++i) {
            const string& contig_id = valid_contigs[i];
            string sequence = assembly.get_sequence(contig_id);
            
            local_processed++;
            ContigRegion region = process_contig(contig_id, sequence, query_kmers, k, threshold_percent);
            
            // only keep regions that meet the threshold
            if (region.kmer_count > 0) {
                local_results.push_back(region);
            }
            
            // progress reporting (only from thread 0)
            if (omp_get_thread_num() == 0 && local_processed % 1000 == 0) {
                uint32_t estimated_total = local_processed * omp_get_num_threads();
                cerr << "progress: ~" << estimated_total << " contigs processed..." << endl;
            }
        }
        
        // merge local results into global results
        #pragma omp critical
        {
            results.insert(results.end(), local_results.begin(), local_results.end());
        }
    }
#else
    // fallback: sequential processing
    uint32_t processed_contigs = 0;
    
    for (const string& contig_id : valid_contigs) {
        string sequence = assembly.get_sequence(contig_id);
        
        processed_contigs++;
        ContigRegion region = process_contig(contig_id, sequence, query_kmers, k, threshold_percent);
        
        // only keep regions that meet the threshold
        if (region.kmer_count > 0) {
            results.push_back(region);
        }
        
        if (processed_contigs % 1000 == 0) {
            cerr << "processed " << processed_contigs << " contigs..." << endl;
        }
    }
#endif
    
    cerr << "processed " << valid_contigs.size() << " contigs, found " 
         << results.size() << " qualifying regions" << endl;
    
    // sort by coverage percentage (descending)
    sort(results.begin(), results.end(), 
         [](const ContigRegion& a, const ContigRegion& b) {
             return a.coverage_percent > b.coverage_percent;
         });
    
    if (!results.empty()) {
        cerr << "best match: contig " << results[0].contig 
             << " (" << results[0].coverage_percent << "% coverage, " 
             << results[0].length << " bp)" << endl;
    }
    
    cerr << "homolog search completed: " << results.size() << " regions found" << endl;
    return results;
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


