#ifndef SEQUENCES_H
#define SEQUENCES_H

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

enum class ResolveSeams {
    NO,              // place Ns for all cases
    READS_ONLY,      // resolve only read seams (requires read sequence data)
    REFERENCE_ONLY,  // resolve only reference seams (requires reference sequence data)
    COMPLETE         // requires both
};

// Convert string to ResolveSeams enum
ResolveSeams string_to_resolve_seams(const std::string& str);

// Convert ResolveSeams enum to string
std::string resolve_seams_to_string(ResolveSeams mode);

// Forward declaration
struct ReadEvent;

// Struct to represent a seam interval that needs sequence resolution
struct SeamInterval {
    std::string lib_id;
    std::string read_id;
    uint32_t start;
    uint32_t end;
    bool reverse_complement;
    bool is_gap;  // true for gaps, false for overlaps
    
    SeamInterval(const std::string& lib_id = "", const std::string& read_id = "",
                uint32_t start = 0, uint32_t end = 0, bool reverse_complement = false, bool is_gap = false)
        : lib_id(lib_id), read_id(read_id), start(start), end(end), 
          reverse_complement(reverse_complement), is_gap(is_gap) {}
};

// Class to hold assembly/reference sequences for all contigs
class AssemblySequences {
private:
    std::unordered_map<std::string, std::string> sequences;
    bool loaded;

public:
    AssemblySequences();
    
    // load sequences from FASTA file
    void load_from_file(const std::string& filename);
    
    // load sequences from map (for R interface)
    void load_from_map(const std::unordered_map<std::string, std::string>& sequence_map);
    
    // get sequence by contig ID
    std::string get_sequence(const std::string& contig_id) const;
    
    // check if sequences are loaded
    bool is_loaded() const { return loaded; }
    
    // get number of sequences
    size_t size() const { return sequences.size(); }
    
    // check if contig exists
    bool has_contig(const std::string& contig_id) const;
    
    // clear all sequences
    void clear();
    
    // get all contig IDs
    std::unordered_set<std::string> get_contig_ids() const;
};

// Class to resolve read intervals from files without loading into memory
class ReadSequences {
private:
    std::unordered_map<std::string, std::string> lib_id_to_filename;  // lib_id -> file path

public:
    ReadSequences();
    
    // register a library file (doesn't load into memory)
    void register_lib_file(const std::string& lib_id, const std::string& filename);
    
    // get sequences for specific read IDs (result passed by reference to avoid copy)
    void get_sequences(const std::string& lib_id, const std::vector<std::string>& read_ids,
                      std::unordered_map<std::string, std::string>& result) const;
    
    // get sequences with caching support
    void get_sequences_with_cache(const std::string& lib_id, const std::vector<std::string>& read_ids,
                                 std::unordered_map<std::string, std::string>& result, const std::string& output_prefix) const;
    
    // utility functions
    void clear();
    bool has_lib(const std::string& lib_id) const;
    std::unordered_set<std::string> get_lib_ids() const;
    
private:
    // caching helper methods
    bool load_sequences_from_cache(const std::string& cache_filename, const std::vector<std::string>& read_ids, 
                                  std::unordered_map<std::string, std::string>& sequences) const;
    void save_sequences_to_cache(const std::string& cache_filename, const std::unordered_map<std::string, std::string>& sequences) const;
};

#endif // SEQUENCES_H
