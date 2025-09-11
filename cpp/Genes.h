#ifndef GENES_H
#define GENES_H

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <unordered_map>

// Forward declarations for Rcpp types (only available when compiling R interface)
#ifdef RCPP_VERSION
#include <Rcpp.h>
#endif

// Structure to represent a single gene
struct Gene {
    std::string gene_id;
    std::string contig;
    uint32_t start;  // 1-based
    uint32_t end;    // 1-based
    char strand;     // '+' or '-'
    std::string desc;
    
    // constructor
    Gene(const std::string& id, const std::string& contig, uint32_t start, 
         uint32_t end, char strand, const std::string& desc = "")
        : gene_id(id), contig(contig), start(start), end(end), strand(strand), desc(desc) {}
};

// Structure for codon translation
struct CodonTable {
    std::string aa_sequence;      // amino acid sequence (64 chars)
    std::string start_sequence;   // start codon markers (64 chars)
    std::string base1_sequence;   // first base of each codon (64 chars)
    std::string base2_sequence;   // second base of each codon (64 chars)  
    std::string base3_sequence;   // third base of each codon (64 chars)
    
    // translate codon to amino acid
    char translate_codon(const std::string& codon) const;
    
    // load from file
    bool load_from_file(const std::string& filepath);
};

// Structure for gene annotation results
struct GeneAnnotation {
    // basic classification
    std::string loc;  // "genic" or "intergenic"
    
    // genic annotations
    std::string gene_id;
    std::string gene_desc;
    uint32_t aa_coord;        // amino acid position in gene (1-based, NA if intergenic)
    std::string variant_codon;
    std::string ref_codon;
    std::string variant_type; // "syn", "non-syn", "ins", "del", "clip"
    std::string mutation_desc; // amino acid change description, e.g., "S:F" for Ser->Phe
    
    // intergenic annotations
    std::string gene_left;
    std::string gene_right;
    std::string orientation_left;   // "upstream" or "downstream"
    std::string orientation_right;  // "upstream" or "downstream"
    uint32_t distance_left;
    uint32_t distance_right;
    
    // constructor for genic
    GeneAnnotation(const std::string& gene_id, const std::string& gene_desc,
                   uint32_t aa_coord, const std::string& variant_codon,
                   const std::string& ref_codon, const std::string& variant_type,
                   const std::string& mutation_desc)
        : loc("genic"), gene_id(gene_id), gene_desc(gene_desc), aa_coord(aa_coord),
          variant_codon(variant_codon), ref_codon(ref_codon), variant_type(variant_type),
          mutation_desc(mutation_desc), distance_left(0), distance_right(0) {}
    
    // constructor for intergenic
    GeneAnnotation(const std::string& gene_left, const std::string& gene_right,
                   const std::string& orientation_left, const std::string& orientation_right,
                   uint32_t distance_left, uint32_t distance_right)
        : loc("intergenic"), aa_coord(0), mutation_desc(""),
          gene_left(gene_left), gene_right(gene_right),
          orientation_left(orientation_left), orientation_right(orientation_right),
          distance_left(distance_left), distance_right(distance_right) {}
};

class Genes {
private:
    // gene storage - flat vector for efficient binary search
    std::vector<Gene> genes;
    
    // gene storage organized by contig for lookup
    std::map<std::string, std::vector<size_t>> contig_to_genes;
    
    // binary search result cache: contig -> position -> gene_index (SIZE_MAX if not genic)
    mutable std::map<std::string, std::map<uint32_t, size_t>> position_search_cache;
    
    // codon table for translation
    CodonTable codon_table;
    
    // codon table sequences (for loading from file)
    std::string aa_sequence;
    std::string base1_sequence;
    std::string base2_sequence;
    std::string base3_sequence;
    
    // reference sequences organized by contig
    std::unordered_map<std::string, std::string> reference_sequences;
    
    bool genes_loaded;
    bool codon_table_loaded;
    bool reference_loaded;
    
    // helper functions
    void sort_genes_by_start();
    void build_contig_to_genes_map();
    size_t binary_search_gene_with_cache(const std::string& contig, uint32_t position) const;
    size_t binary_search_gene(const std::string& contig, uint32_t position) const;
    GeneAnnotation find_flanking_genes(const std::string& contig, uint32_t position) const;
    uint32_t get_aa_coordinate(uint32_t position, uint32_t gene_start, char strand, 
                              uint32_t gene_end) const;
    
    // reference sequence helper functions
    std::string get_reference_codon(const std::string& contig, uint32_t position, 
                                   const Gene& gene) const;
    std::string apply_mutation_to_codon(const std::string& ref_codon, int codon_position, 
                                       char variant_base) const;
    std::string reverse_complement(const std::string& sequence) const;
    
public:
    Genes();
    
    // load gene table from file
    bool load_gene_table(const std::string& filepath);
    
    // load codon table from file or standard name
    bool load_codon_table(const std::string& filepath_or_name);
    
    // load reference sequences from FASTA file
    bool load_reference_sequences(const std::string& fasta_path);
    
    // load gene table from C++ vector
    bool load_gene_table_from_vector(const std::vector<Gene>& gene_vector);
    
    // load reference sequences from C++ map
    bool load_reference_sequences_from_map(const std::unordered_map<std::string, std::string>& sequences_map);
    
    // main annotation function
    GeneAnnotation annotate_variant(const std::string& contig, uint32_t position,
                                   const std::string& mutation_type, 
                                   const std::string& mutation_sequence) const;
    
    // check if genes are loaded
    bool is_loaded() const { 
        return genes_loaded && codon_table_loaded && reference_loaded;
    }
    
    // check if reference sequences are loaded
    bool has_reference_sequences() const { return reference_loaded; }
    
    // get gene count
    size_t get_gene_count() const;
    size_t get_gene_count(const std::string& contig) const;
};

#endif // GENES_H
