#ifndef HOMOLOGS_H
#define HOMOLOGS_H

#include <string>
#include <vector>
#include <unordered_map>

struct ContigRegion {
    std::string contig;
    uint32_t start;
    uint32_t end;
    uint32_t length;
    uint32_t kmer_count;
    double coverage_percent;
};

class Homologs {
private:
    int num_threads;
    
    // extract kmers evenly spaced from query region
    std::vector<std::string> extract_query_kmers(const std::string& sequence, uint32_t start, uint32_t end, 
                                                 uint32_t k, uint32_t num_kmers);
    
    // build hash table of kmers for a single contig
    std::unordered_map<std::string, std::vector<uint32_t>> build_contig_kmer_index(const std::string& sequence, uint32_t k);
    
    // process a single contig for homologous regions
    ContigRegion process_contig(const std::string& contig_id,
                               const std::string& sequence,
                               const std::vector<std::string>& query_kmers,
                               uint32_t k,
                               double threshold_percent);

public:
    Homologs(int num_threads = 0);
    ~Homologs();
    
    // main search function (file-based)
    std::vector<ContigRegion> search_homologs(const std::string& fasta_file,
                                             const std::string& query_contig,
                                             uint32_t query_start,
                                             uint32_t query_end,
                                             uint32_t k,
                                             uint32_t num_kmers,
                                             double threshold_percent);
    
    // main search function (pre-loaded sequences)
    std::vector<ContigRegion> search_homologs(const class AssemblySequences& assembly,
                                             const std::string& query_contig,
                                             uint32_t query_start,
                                             uint32_t query_end,
                                             uint32_t k,
                                             uint32_t num_kmers,
                                             double threshold_percent);
    
    // write results to file
    void write_results(const std::vector<ContigRegion>& regions, const std::string& output_file);
};

#endif // HOMOLOGS_H
