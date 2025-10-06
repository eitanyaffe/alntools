#include "Genes.h"
#include "utils.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cctype>
#include <unordered_set>

using namespace std;

// CodonTable implementation
char CodonTable::translate_codon(const std::string& codon) const {
    if (codon.length() != 3) return 'X'; // invalid codon
    
    // find codon index in the standard genetic code order
    int index = -1;
    for (int i = 0; i < 64; i++) {
        if (base1_sequence[i] == toupper(codon[0]) &&
            base2_sequence[i] == toupper(codon[1]) &&
            base3_sequence[i] == toupper(codon[2])) {
            index = i;
            break;
        }
    }
    
    if (index >= 0 && index < 64) {
        return aa_sequence[index];
    }
    return 'X'; // unknown
}

bool CodonTable::load_from_file(const std::string& filepath) {
    ifstream file(filepath);
    if (!file.is_open()) {
        cerr << "error: cannot open codon table file: " << filepath << endl;
        return false;
    }
    
    string line;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        size_t eq_pos = line.find('=');
        if (eq_pos == string::npos) continue;
        
        string key = line.substr(0, eq_pos);
        string value = line.substr(eq_pos + 1);
        
        // trim whitespace from key
        size_t key_start = key.find_first_not_of(" \t\r\n");
        size_t key_end = key.find_last_not_of(" \t\r\n");
        if (key_start != string::npos && key_end != string::npos) {
            key = key.substr(key_start, key_end - key_start + 1);
        }
        
        // trim whitespace from value
        size_t start = value.find_first_not_of(" \t\r\n");
        size_t end = value.find_last_not_of(" \t\r\n");
        if (start != string::npos && end != string::npos) {
            value = value.substr(start, end - start + 1);
        }
        
        if (key == "AA") {
            aa_sequence = value;
        } else if (key == "START") {
            start_sequence = value;
        } else if (key == "Base1") {
            base1_sequence = value;
        } else if (key == "Base2") {
            base2_sequence = value;
        } else if (key == "Base3") {
            base3_sequence = value;
        }
    }
    
    // validate that all sequences are loaded and have correct length
    if (aa_sequence.length() != 64 || start_sequence.length() != 64 ||
        base1_sequence.length() != 64 || base2_sequence.length() != 64 ||
        base3_sequence.length() != 64) {
        cerr << "error: invalid codon table format in " << filepath << endl;
        return false;
    }
    
    cout << "loaded codon table from " << filepath << endl;
    return true;
}

// Genes implementation
Genes::Genes() : genes_loaded(false), codon_table_loaded(false), reference_loaded(false) {
}

bool Genes::load_gene_table(const std::string& filepath) {
    cout << "loading gene table from " << filepath << endl;
    
    ifstream file(filepath);
    if (!file.is_open()) {
        cerr << "error: cannot open gene table file: " << filepath << endl;
        return false;
    }
    
    string line;
    bool header_read = false;
    
    while (getline(file, line)) {
        if (line.empty()) continue;
        
        // skip header
        if (!header_read) {
            header_read = true;
            continue;
        }
        
        istringstream ss(line);
        string gene_id, contig, strand_str, desc;
        uint32_t start, end;
        
        // parse tab-delimited: gene, contig, start, end, strand, [other fields], desc
        if (!(ss >> gene_id >> contig >> start >> end >> strand_str)) {
            cerr << "warning: skipping malformed line: " << line << endl;
            continue;
        }
        
        char strand = '+';
        if (!strand_str.empty()) {
            strand = strand_str[0];
        }
        
        // read the rest of the line as description (after the 5 required fields)
        string remaining_line;
        getline(ss, remaining_line);
        
        // trim leading whitespace and use as description
        size_t desc_start = remaining_line.find_first_not_of(" \t");
        if (desc_start != string::npos) {
            desc = remaining_line.substr(desc_start);
        }
        
        // add gene
        genes.emplace_back(gene_id, contig, start, end, strand, desc);
    }
    
    // finalize gene loading with shared code
    finalize_gene_loading();
    return true;
}

bool Genes::load_gene_table_from_vector(const std::vector<Gene>& gene_vector) {
    cout << "loading gene table from vector" << endl;
    
    // clear existing genes
    genes.clear();
    contig_to_genes.clear();
    position_search_cache.clear();
    genes_loaded = false;
    
    // copy genes from vector
    genes = gene_vector;
    
    cout << "processing " << genes.size() << " genes from vector" << endl;
    
    // finalize gene loading with shared code
    finalize_gene_loading();
    return true;
}

bool Genes::load_codon_table(const std::string& filepath_or_name) {
    string filepath = filepath_or_name;
    
    // if it doesn't contain a path separator, assume it's a standard table name
    if (filepath.find('/') == string::npos && filepath.find('\\') == string::npos) {
        // construct path relative to codon_tables directory
        filepath = "codon_tables/" + filepath_or_name;
    }
    
    codon_table_loaded = codon_table.load_from_file(filepath);
    return codon_table_loaded;
}


// sort_genes_by_start implementation moved to end of file

size_t Genes::binary_search_gene(const std::string& contig, uint32_t position) const {
    auto contig_it = contig_to_genes.find(contig);
    if (contig_it == contig_to_genes.end()) {
        return SIZE_MAX; // contig not found
    }
    
    const vector<size_t>& gene_indices = contig_it->second;
    
    // binary search for gene containing position
    size_t left = 0;
    size_t right = gene_indices.size();
    
    while (left < right) {
        size_t mid = left + (right - left) / 2;
        const Gene& gene = genes[gene_indices[mid]];
        
        if (gene.end < position) {
            left = mid + 1;
        } else if (gene.start > position) {
            right = mid;
        } else {
            // position is within gene[mid]
            return gene_indices[mid];
        }
    }
    
    return SIZE_MAX; // not found
}

size_t Genes::binary_search_gene_with_cache(const std::string& contig, uint32_t position) const {
    // check cache first
    auto contig_cache_it = position_search_cache.find(contig);
    if (contig_cache_it != position_search_cache.end()) {
        auto pos_it = contig_cache_it->second.find(position);
        if (pos_it != contig_cache_it->second.end()) {
            return pos_it->second;
        }
    }
    
    // cache miss - do binary search
    size_t result = binary_search_gene(contig, position);
    
    // store result in cache
    position_search_cache[contig][position] = result;
    
    return result;
}


uint32_t Genes::get_aa_coordinate(uint32_t position, uint32_t gene_start, char strand, uint32_t gene_end) const {
    if (strand == '+') {
        return ((position - gene_start) / 3) + 1;
    } else {
        return ((gene_end - position) / 3) + 1;
    }
}

GeneAnnotation Genes::annotate_variant(const std::string& contig, uint32_t position,
                                      const std::string& mutation_type, 
                                      const std::string& mutation_sequence) const {
    massert(is_loaded(), "Genes::annotate_variant called before genes loaded");
    
    // use cached binary search to find gene
    size_t gene_idx = binary_search_gene_with_cache(contig, position);
    
    if (gene_idx != SIZE_MAX) {
        // position is genic
        const Gene& gene = genes[gene_idx];
        
        // calculate amino acid coordinate
        uint32_t aa_coord = get_aa_coordinate(position, gene.start, gene.strand, gene.end);
        
        // determine variant type and mutation description using mutation data directly
        string var_type = mutation_type;
        string mutation_desc = "";
        string ref_codon = "NNN";  // placeholder
        string var_codon = "NNN";  // placeholder
        
        if (mutation_type == "sub" && mutation_sequence.length() >= 2) {
            // For substitutions: mutation_sequence[0] = variant base, mutation_sequence[1] = reference base
            char variant_base = mutation_sequence[0];
            char reference_base = mutation_sequence[1];
            
            // extract proper codon context from reference sequence
            ref_codon = get_reference_codon(contig, position, gene);
            
            if (ref_codon != "NNN") {
                // calculate position within codon (0, 1, or 2)
                int codon_position;
                if (gene.strand == '+') {
                    codon_position = (position - gene.start) % 3;
                } else {
                    codon_position = (gene.end - position) % 3;
                }
                
                // apply mutation to get variant codon
                char mutation_base = (gene.strand == '+') ? variant_base : 
                                    (variant_base == 'A' ? 'T' : variant_base == 'T' ? 'A' : 
                                     variant_base == 'G' ? 'C' : variant_base == 'C' ? 'G' : 'N');
                var_codon = apply_mutation_to_codon(ref_codon, codon_position, mutation_base);
                
                // translate codons to amino acids
                char ref_aa = codon_table.translate_codon(ref_codon);
                char var_aa = codon_table.translate_codon(var_codon);
                
                var_type = (ref_aa == var_aa) ? "syn" : "non-syn";
                
                // generate mutation description
                if (ref_aa != 'X' && var_aa != 'X') {
                    mutation_desc = string(1, ref_aa) + ":" + string(1, var_aa);
                } else {
                    mutation_desc = string(1, reference_base) + ":" + string(1, variant_base);
                }
            } else {
                // codon extraction failed - this should not happen with proper reference
                cerr << "error: failed to extract codon for position " << position << " on " << contig << endl;
                ref_codon = "NNN";
                var_codon = "NNN";
                mutation_desc = string(1, reference_base) + ":" + string(1, variant_base);
            }
        } else if (mutation_type == "ins") {
            mutation_desc = "ins";
            ref_codon = "";
            var_codon = mutation_sequence;  // inserted bases
        } else if (mutation_type == "del") {
            mutation_desc = "del";
            ref_codon = mutation_sequence;  // deleted bases
            var_codon = "";
        } else {
            mutation_desc = mutation_type; // fallback for other types
        }
        
        return GeneAnnotation(gene.gene_id, gene.desc, aa_coord, var_codon, ref_codon, var_type, mutation_desc);
    }
    
    // position is intergenic - find adjacent genes
    auto contig_genes_it = contig_to_genes.find(contig);
    if (contig_genes_it == contig_to_genes.end()) {
        return GeneAnnotation("", "", "", "", 0, 0); // no genes on this contig
    }
    
    const auto& gene_indices = contig_genes_it->second;
    
    // find left and right adjacent genes
    string gene_left = "", gene_right = "";
    string orient_left = "", orient_right = "";
    uint32_t dist_left = 0, dist_right = 0;
    
    // find closest gene to the left
    for (int i = gene_indices.size() - 1; i >= 0; i--) {
        const Gene& gene = genes[gene_indices[i]];
        if (gene.end < position) {
            gene_left = gene.gene_id;
            dist_left = position - gene.end;
            orient_left = (gene.strand == '+') ? "downstream" : "upstream";
            break;
        }
    }
    
    // find closest gene to the right
    for (size_t i = 0; i < gene_indices.size(); i++) {
        const Gene& gene = genes[gene_indices[i]];
        if (gene.start > position) {
            gene_right = gene.gene_id;
            dist_right = gene.start - position;
            orient_right = (gene.strand == '+') ? "upstream" : "downstream";
            break;
        }
    }
    
    return GeneAnnotation(gene_left, gene_right, orient_left, orient_right, dist_left, dist_right);
}

size_t Genes::get_gene_count() const {
    return genes.size();
}

size_t Genes::get_gene_count(const std::string& contig) const {
    auto it = contig_to_genes.find(contig);
    return (it != contig_to_genes.end()) ? it->second.size() : 0;
}

bool Genes::load_reference_sequences(const std::string& fasta_path) {
    cout << "loading reference sequences from " << fasta_path << endl;
    
    // clear any existing sequences
    reference_sequences.clear();
    reference_loaded = false;
    
    try {
        // use empty set to load all contigs
        unordered_set<string> empty_set;
        read_fasta(fasta_path, empty_set, reference_sequences);
        
        reference_loaded = true;
        cout << "loaded " << reference_sequences.size() << " reference sequences" << endl;
        return true;
    } catch (const exception& e) {
        cerr << "error loading reference sequences: " << e.what() << endl;
        return false;
    }
}

bool Genes::load_reference_sequences_from_map(const std::unordered_map<std::string, std::string>& sequences_map) {
    cout << "loading reference sequences from map" << endl;
    
    // clear any existing sequences
    reference_sequences.clear();
    reference_loaded = false;
    
    // copy sequences from map
    reference_sequences = sequences_map;
    
    reference_loaded = true;
    cout << "loaded " << reference_sequences.size() << " reference sequences" << endl;
    return true;
}

std::string Genes::reverse_complement(const std::string& sequence) const {
    string result;
    result.reserve(sequence.length());
    
    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        switch (toupper(*it)) {
            case 'A': result += 'T'; break;
            case 'T': result += 'A'; break;
            case 'G': result += 'C'; break;
            case 'C': result += 'G'; break;
            default: result += 'N'; break;  // handle ambiguous bases
        }
    }
    return result;
}

std::string Genes::get_reference_codon(const std::string& contig, uint32_t position, 
                                      const Gene& gene) const {
    if (!reference_loaded) {
        cerr << "error: reference sequences not loaded - cannot extract codon" << endl;
        return "NNN";
    }
    
    auto seq_it = reference_sequences.find(contig);
    if (seq_it == reference_sequences.end()) {
        cerr << "error: contig " << contig << " not found in reference sequences" << endl;
        return "NNN";
    }
    
    const string& sequence = seq_it->second;
    
    // calculate codon start position (0-based for string indexing)
    uint32_t codon_start;
    if (gene.strand == '+') {
        // forward strand: find start of codon containing this position
        uint32_t offset_from_gene_start = position - gene.start;
        codon_start = gene.start + (offset_from_gene_start / 3) * 3 - 1;  // convert to 0-based
    } else {
        // reverse strand: find start of codon containing this position
        uint32_t offset_from_gene_end = gene.end - position;
        codon_start = gene.end - (offset_from_gene_end / 3) * 3 - 3;  // convert to 0-based, go back 3
    }
    
    // extract 3 bases
    if (codon_start + 2 >= sequence.length()) {
        return "NNN";  // out of bounds
    }
    
    string codon = sequence.substr(codon_start, 3);
    
    // reverse complement if on minus strand
    if (gene.strand == '-') {
        codon = reverse_complement(codon);
    }
    
    return codon;
}

std::string Genes::apply_mutation_to_codon(const std::string& ref_codon, int codon_position, 
                                          char variant_base) const {
    if (ref_codon.length() != 3 || codon_position < 0 || codon_position > 2) {
        return ref_codon;  // invalid input, return unchanged
    }
    
    string variant_codon = ref_codon;
    variant_codon[codon_position] = toupper(variant_base);
    return variant_codon;
}

void Genes::sort_genes_by_start() {
    sort(genes.begin(), genes.end(), [](const Gene& a, const Gene& b) {
        if (a.contig != b.contig) {
            return a.contig < b.contig;
        }
        return a.start < b.start;
    });
}

void Genes::build_contig_to_genes_map() {
    contig_to_genes.clear();
    
    for (size_t i = 0; i < genes.size(); i++) {
        const Gene& gene = genes[i];
        contig_to_genes[gene.contig].push_back(i);
    }
}

void Genes::finalize_gene_loading() {
    // count unique contigs
    unordered_set<string> contigs_with_genes;
    for (const auto& gene : genes) {
        contigs_with_genes.insert(gene.contig);
    }
    
    // sort genes by start position for binary search
    sort_genes_by_start();
    
    // build contig to gene index map
    build_contig_to_genes_map();
    
    genes_loaded = true;
    cout << "loaded " << genes.size() << " genes from " << contigs_with_genes.size() << " contigs" << endl;
}
