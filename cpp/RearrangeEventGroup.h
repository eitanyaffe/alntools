#ifndef REARRANGE_EVENT_GROUP_H
#define REARRANGE_EVENT_GROUP_H

#include "RearrangeReadEvent.h"
#include <vector>
#include <map>
#include <string>

// representative event (no read_strand field)
struct Event {
    std::string event_id;
    RearrangementType type;
    std::string contig_id;
    uint32_t out_clip;
    uint32_t in_clip;
    std::string element_contig;
    std::string element_strand;
    uint32_t element_start;
    uint32_t element_end;
    std::string read_seams;      // seams in read coordinates and sequences
    std::string assembly_seams;  // seams in assembly coordinates and sequences
    
    Event(const std::string& event_id = "", RearrangementType type = RearrangementType::LARGE_DELETE,
          const std::string& contig_id = "",
          uint32_t out_clip = 0, uint32_t in_clip = 0,
          const std::string& element_contig = "", const std::string& element_strand = "",
          uint32_t element_start = 0, uint32_t element_end = 0,
          const std::string& read_seams = "", const std::string& assembly_seams = "")
        : event_id(event_id), type(type), contig_id(contig_id),
          out_clip(out_clip), in_clip(in_clip),
          element_contig(element_contig), element_strand(element_strand),
          element_start(element_start), element_end(element_end),
          read_seams(read_seams), assembly_seams(assembly_seams) {}
    
    // create key for exact matching (all fields)
    std::string create_key() const;
};

// simple Union-Find data structure for connected components
class UnionFind {
private:
    std::vector<int> parent;
    std::vector<int> rank;

public:
    UnionFind(int n);
    
    int find(int x);
    void unite(int x, int y);
    std::vector<std::vector<int>> get_components();
};

// worker class for grouping similar events across libraries
class RearrangeEventGroup {
private:
    int max_margin;
    
    // intermediate: flat list of all read events with library info
    struct IndexedReadEvent {
        std::string lib_id;
        size_t event_index;
        const ReadEvent* event;
        
        IndexedReadEvent(const std::string& lib_id, size_t event_index, const ReadEvent* event)
            : lib_id(lib_id), event_index(event_index), event(event) {}
    };

public:
    RearrangeEventGroup(int max_margin = 10);
    
    // main grouping function - pass all containers by reference
    void group_events(const std::map<std::string, std::vector<ReadEvent>>& lib_events,
                     std::vector<Event>& representative_events,
                     std::vector<std::vector<int>>& event_groups,
                     std::map<std::string, std::map<std::string, int>>& lib_event_support);

private:
    // similarity comparison
    bool is_similar(const ReadEvent& a, const ReadEvent& b) const;
    
    // exact match comparison (for representative selection)
    bool is_exact_match(const ReadEvent& a, const ReadEvent& b) const;
    
    // create representative event from read event
    Event create_representative_event(const ReadEvent& read_event, const std::string& event_id) const;
    
    // flatten library events into indexed list
    void flatten_events(const std::map<std::string, std::vector<ReadEvent>>& lib_events,
                       std::vector<IndexedReadEvent>& all_events);
    
    // build similarity graph and find connected components
    std::vector<std::vector<int>> find_connected_components(const std::vector<IndexedReadEvent>& all_events);
    
    // select representative event for each group
    void select_representatives(const std::vector<IndexedReadEvent>& all_events,
                               const std::vector<std::vector<int>>& event_groups,
                               std::vector<Event>& representative_events);
    
    // calculate support counts per library
    void calculate_support_counts(const std::map<std::string, std::vector<ReadEvent>>& lib_events,
                                 const std::vector<IndexedReadEvent>& all_events,
                                 const std::vector<std::vector<int>>& event_groups,
                                 const std::vector<Event>& representative_events,
                                 std::map<std::string, std::map<std::string, int>>& lib_event_support);
};

#endif // REARRANGE_EVENT_GROUP_H
