#ifndef HOMOLOGS_H
#define HOMOLOGS_H

#include <string>
#include <vector>
#include <unordered_map>

struct KmerHit {
    std::string contig;
    uint32_t position;
};

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
    // extract kmers evenly spaced from query region
    std::vector<std::string> extract_query_kmers(const std::string& sequence, uint32_t start, uint32_t end, 
                                                 uint32_t k, uint32_t num_kmers);
    
    // build hash table of all kmers in assembly
    std::unordered_map<std::string, std::vector<KmerHit>> build_kmer_index(const class AssemblySequences& assembly, uint32_t k);
    
    // find regions with sufficient kmer coverage
    std::vector<ContigRegion> find_homologous_regions(const std::vector<std::string>& query_kmers,
                                                     const std::unordered_map<std::string, std::vector<KmerHit>>& kmer_index,
                                                     double threshold_percent);

public:
    Homologs();
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
