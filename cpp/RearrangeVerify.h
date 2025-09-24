#ifndef REARRANGE_VERIFY_H
#define REARRANGE_VERIFY_H

#include <string>
#include <unordered_map>

// Forward declarations
struct ReadEvent;

class RearrangeVerify {
private:
    std::unordered_map<std::string, std::string> contig_sequences;
    std::unordered_map<std::string, std::string> read_sequences;

public:
    RearrangeVerify() = default;
    
    // load sequences from maps (for R interface)
    void load_sequences_from_maps(const std::unordered_map<std::string, std::string>& contigs,
                                 const std::unordered_map<std::string, std::string>& reads);
    
    // load sequences from files
    void load_contig_sequences_from_file(const std::string& filename);
    void load_read_sequences_from_file(const std::string& filename);
    
    // verify a read event
    bool verify_event(const ReadEvent& event) const;
    
    // clear loaded sequences
    void clear();
    
private:
    // helper functions for sequence validation
    bool validate_anchors(const ReadEvent& event) const;
    bool validate_element(const ReadEvent& event) const;
    bool validate_shims(const ReadEvent& event) const;
};

#endif // REARRANGE_VERIFY_H
