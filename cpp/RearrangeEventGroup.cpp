#include "RearrangeEventGroup.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>

using namespace std;

// Event implementation
string Event::create_key() const
{
    ostringstream oss;
    oss << static_cast<int>(type) << "|"
        << contig_id << "|"
        << out_clip << "|"
        << in_clip << "|"
        << element_contig << "|"
        << element_strand << "|"
        << element_start << "|"
        << element_end << "|"
        << left_shim << "|"
        << right_shim << "|"
        << middle_shim;
    return oss.str();
}

// UnionFind implementation
UnionFind::UnionFind(int n) : parent(n), rank(n, 0)
{
    for (int i = 0; i < n; ++i) {
        parent[i] = i;
    }
}

int UnionFind::find(int x)
{
    if (parent[x] != x) {
        parent[x] = find(parent[x]); // path compression
    }
    return parent[x];
}

void UnionFind::unite(int x, int y)
{
    int px = find(x);
    int py = find(y);
    
    if (px == py) return;
    
    // union by rank
    if (rank[px] < rank[py]) {
        parent[px] = py;
    } else if (rank[px] > rank[py]) {
        parent[py] = px;
    } else {
        parent[py] = px;
        rank[px]++;
    }
}

vector<vector<int>> UnionFind::get_components()
{
    map<int, vector<int>> root_to_members;
    
    for (int i = 0; i < static_cast<int>(parent.size()); ++i) {
        int root = find(i);
        root_to_members[root].push_back(i);
    }
    
    vector<vector<int>> components;
    for (const auto& pair : root_to_members) {
        components.push_back(pair.second);
    }
    
    return components;
}

// RearrangeEventGroup implementation
RearrangeEventGroup::RearrangeEventGroup(int max_margin) : max_margin(max_margin) {}

void RearrangeEventGroup::group_events(const map<string, vector<ReadEvent>>& lib_events,
                                       vector<Event>& representative_events,
                                       vector<vector<int>>& event_groups,
                                       map<string, map<string, int>>& lib_event_support)
{
    cout << "grouping events across libraries..." << endl;
    
    // clear output containers
    representative_events.clear();
    event_groups.clear();
    lib_event_support.clear();
    
    // step 1: flatten all events into indexed list
    vector<IndexedReadEvent> all_events;
    flatten_events(lib_events, all_events);
    
    if (all_events.empty()) {
        cout << "no events to group" << endl;
        return;
    }
    
    cout << "total events across all libraries: " << all_events.size() << endl;
    
    // step 2: find connected components based on similarity
    event_groups = find_connected_components(all_events);
    
    cout << "found " << event_groups.size() << " event groups" << endl;
    
    // step 3: select representative event for each group
    select_representatives(all_events, event_groups, representative_events);
    
    // step 4: calculate support counts per library
    calculate_support_counts(lib_events, all_events, event_groups, representative_events, lib_event_support);
    
    cout << "event grouping completed" << endl;
}


bool RearrangeEventGroup::is_similar(const ReadEvent& a, const ReadEvent& b) const
{
    // must be same type
    if (a.type != b.type) return false;
    
    // must be same contig
    if (a.contig_id != b.contig_id) return false;
    
    // coordinates must be within margin
    if (abs(static_cast<int>(a.out_clip) - static_cast<int>(b.out_clip)) > max_margin) return false;
    if (abs(static_cast<int>(a.in_clip) - static_cast<int>(b.in_clip)) > max_margin) return false;
    
    // element details must match (if present)
    if (a.element_contig != b.element_contig) return false;
    if (a.element_strand != b.element_strand) return false;
    
    if (a.element_start != 0 && b.element_start != 0) {
        if (abs(static_cast<int>(a.element_start) - static_cast<int>(b.element_start)) > max_margin) return false;
        if (abs(static_cast<int>(a.element_end) - static_cast<int>(b.element_end)) > max_margin) return false;
    }
    
    // shim lengths must be similar (within margin)
    if (abs(static_cast<int>(a.left_shim.length()) - static_cast<int>(b.left_shim.length())) > max_margin) return false;
    if (abs(static_cast<int>(a.right_shim.length()) - static_cast<int>(b.right_shim.length())) > max_margin) return false;
    if (abs(static_cast<int>(a.middle_shim.length()) - static_cast<int>(b.middle_shim.length())) > max_margin) return false;
    
    // note: we ignore read_strand as requested
    // note: we ignore shim sequences, only check lengths
    
    return true;
}

bool RearrangeEventGroup::is_exact_match(const ReadEvent& a, const ReadEvent& b) const
{
    // all key fields must be identical (except read_strand)
    return a.type == b.type &&
           a.contig_id == b.contig_id &&
           a.out_clip == b.out_clip &&
           a.in_clip == b.in_clip &&
           a.element_contig == b.element_contig &&
           a.element_strand == b.element_strand &&
           a.element_start == b.element_start &&
           a.element_end == b.element_end &&
           a.left_shim == b.left_shim &&
           a.right_shim == b.right_shim &&
           a.middle_shim == b.middle_shim;
}

Event RearrangeEventGroup::create_representative_event(const ReadEvent& read_event, const string& event_id) const
{
    return Event(event_id, read_event.type, read_event.contig_id,
                read_event.out_clip, read_event.in_clip,
                read_event.element_contig, read_event.element_strand,
                read_event.element_start, read_event.element_end,
                read_event.left_shim, read_event.right_shim, read_event.middle_shim);
}

void RearrangeEventGroup::flatten_events(const map<string, vector<ReadEvent>>& lib_events,
                                         vector<IndexedReadEvent>& all_events)
{
    all_events.clear();
    // print the number of libraries in lib_events
    cout << "flatten_events, number of events: " << lib_events.size() << " libraries" << endl;
    for (const auto& lib_pair : lib_events) {
        const string& lib_id = lib_pair.first;
        const vector<ReadEvent>& events = lib_pair.second;
        
        for (size_t i = 0; i < events.size(); ++i) {
            all_events.emplace_back(lib_id, i, &events[i]);
        }
    }
}

vector<vector<int>> RearrangeEventGroup::find_connected_components(const vector<IndexedReadEvent>& all_events)
{
    const int n = static_cast<int>(all_events.size());
    UnionFind uf(n);
    
    // build similarity graph by comparing all pairs
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (is_similar(*all_events[i].event, *all_events[j].event)) {
                uf.unite(i, j);
            }
        }
    }
    
    return uf.get_components();
}

void RearrangeEventGroup::select_representatives(const vector<IndexedReadEvent>& all_events,
                                                const vector<vector<int>>& event_groups,
                                                vector<Event>& representative_events)
{
    representative_events.clear();
    representative_events.reserve(event_groups.size());
    
    // create representatives directly for each group
    for (const vector<int>& group : event_groups) {
        // count exact matches for each unique event in group
        map<string, pair<size_t, int>> event_key_to_info; // key -> (first_index, count)
        
        for (int event_idx : group) {
            const ReadEvent* event = all_events[event_idx].event;
            
            // create temporary representative to get key
            Event temp_event = create_representative_event(*event, "temp");
            string key = temp_event.create_key();
            
            if (event_key_to_info.find(key) == event_key_to_info.end()) {
                event_key_to_info[key] = make_pair(event_idx, 1);
            } else {
                event_key_to_info[key].second++;
            }
        }
        
        // find event with highest exact match count
        size_t best_event_idx = group[0];
        int best_count = 0;
        
        for (const auto& pair : event_key_to_info) {
            if (pair.second.second > best_count) {
                best_count = pair.second.second;
                best_event_idx = pair.second.first;
            }
        }
        
        // create representative event without ID
        const ReadEvent* best_event = all_events[best_event_idx].event;
        Event rep_event = create_representative_event(*best_event, "temp");
        representative_events.push_back(rep_event);
    }
    
    // create temporary container with indices for sorting to determine ID assignment
    vector<pair<size_t, Event*>> temp_for_sorting;
    for (size_t i = 0; i < representative_events.size(); ++i) {
        temp_for_sorting.emplace_back(i, &representative_events[i]);
    }
    
    // sort by event properties to determine ID order
    sort(temp_for_sorting.begin(), temp_for_sorting.end(), 
         [](const pair<size_t, Event*>& a, const pair<size_t, Event*>& b) {
             if (a.second->contig_id != b.second->contig_id) {
                 return a.second->contig_id < b.second->contig_id;
             }
             if (a.second->out_clip != b.second->out_clip) {
                 return a.second->out_clip < b.second->out_clip;
             }
             return a.second->in_clip < b.second->in_clip;
         });
    
    // assign sequential IDs based on sorted order
    for (size_t i = 0; i < temp_for_sorting.size(); ++i) {
        temp_for_sorting[i].second->event_id = "ev" + to_string(i + 1);
    }
}

void RearrangeEventGroup::calculate_support_counts(const map<string, vector<ReadEvent>>& lib_events,
                                                   const vector<IndexedReadEvent>& all_events,
                                                   const vector<vector<int>>& event_groups,
                                                   const vector<Event>& representative_events,
                                                   map<string, map<string, int>>& lib_event_support)
{
    lib_event_support.clear();
    
    // initialize support counts
    for (const auto& lib_pair : lib_events) {
        const string& lib_id = lib_pair.first;
        lib_event_support[lib_id] = map<string, int>();
        
        for (const Event& event : representative_events) {
            lib_event_support[lib_id][event.event_id] = 0;
        }
    }
    
    // count support for each event in each library
    for (size_t group_idx = 0; group_idx < event_groups.size(); ++group_idx) {
        const string& event_id = representative_events[group_idx].event_id;
        const vector<int>& group = event_groups[group_idx];
        
        // count events per library in this group
        map<string, int> lib_counts;
        for (int event_idx : group) {
            const string& lib_id = all_events[event_idx].lib_id;
            lib_counts[lib_id]++;
        }
        
        // update support matrix
        for (const auto& count_pair : lib_counts) {
            lib_event_support[count_pair.first][event_id] = count_pair.second;
        }
    }
}
