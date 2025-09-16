#include "RearrangeManager.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

string SetRearrangementData::create_key() const
{
    ostringstream oss;
    oss << static_cast<int>(type) << "|"
        << contig_id << "|"
        << strand << "|"
        << out_clip << "|"
        << in_clip << "|"
        << element_contig << "|"
        << element_strand << "|"
        << element_start << "|"
        << element_end;
    return oss.str();
}

string RearrangementOutputRow::create_key() const
{
    ostringstream oss;
    oss << static_cast<int>(type) << "|"
        << contig_id << "|"
        << strand << "|"
        << out_clip << "|"
        << in_clip << "|"
        << element_contig << "|"
        << element_strand << "|"
        << element_start << "|"
        << element_end;
    return oss.str();
}

// RearrangeManager Implementation
RearrangeManager::RearrangeManager(const std::map<std::string, AlignmentStore>& stores,
                                   const std::vector<Interval>& intervals,
                                   VerifyRearrange* verifier,
                                   int max_gap,
                                   int min_element_length,
                                   int min_anchor_length,
                                   double max_mutations_percent)
    : stores(stores)
    , intervals(intervals)
    , max_gap(max_gap)
    , min_element_length(min_element_length)
    , min_anchor_length(min_anchor_length)
    , max_mutations_percent(max_mutations_percent)
    , verifier(verifier)
{
    // extract ordered library IDs and create rearrangers
    for (const auto& entry : stores) {
        library_ids.push_back(entry.first);
        rearrangers[entry.first] = std::make_unique<Rearrange>(
            entry.second, verifier, max_gap, min_element_length, min_anchor_length, max_mutations_percent);
    }
}

void RearrangeManager::execute()
{
    cout << "executing rearrangement detection across " << library_ids.size() << " libraries..." << endl;
    
    // Pass 1: Collect all event keys
    collect_all_event_keys();
    
    // Sort keys and assign event IDs
    sort_and_assign_event_ids();
    
    // Pass 2: Process libraries with assigned event IDs
    process_libraries_with_event_ids();
    
    // Cross-library aggregation and matrices
    aggregate_across_libraries();
    
    cout << "rearrangement detection completed" << endl;
}

void RearrangeManager::collect_all_event_keys()
{
    cout << "collecting event keys from all libraries..." << endl;
    
    all_event_keys.clear();
    
    for (const std::string& lib_id : library_ids) {
        cout << "processing library: " << lib_id << endl;
        
        // ensure read-to-alignment index is built
        const AlignmentStore& store = stores.at(lib_id);
        if (!store.is_read_alignment_index_built()) {
            cerr << "error: read alignment index not built for library " << lib_id << endl;
            continue;
        }
        
        // set interval filtering for this library
        rearrangers[lib_id]->set_reads_and_alignments(intervals);
        
        // get events for this library
        rearrangers[lib_id]->get_events();
        
        // collect keys from this library
        std::set<std::string> lib_keys = rearrangers[lib_id]->get_event_keys();
        all_event_keys.insert(lib_keys.begin(), lib_keys.end());
        
        cout << "library " << lib_id << ": found " << lib_keys.size() << " unique event types" << endl;
    }
    
    cout << "collected " << all_event_keys.size() << " unique event keys across all libraries" << endl;
}

void RearrangeManager::sort_and_assign_event_ids()
{
    cout << "sorting and assigning event IDs..." << endl;
    
    // convert set to vector for sorting
    std::vector<std::string> sorted_keys(all_event_keys.begin(), all_event_keys.end());
    
    // parse keys and sort by contig, then clip_out
    std::sort(sorted_keys.begin(), sorted_keys.end(), [](const std::string& a, const std::string& b) {
        // keys are: type|contig_id|strand|out_clip|in_clip|element_contig|element_strand|element_start|element_end
        std::vector<std::string> parts_a, parts_b;
        std::stringstream ss_a(a), ss_b(b);
        std::string item;
        
        while (std::getline(ss_a, item, '|')) parts_a.push_back(item);
        while (std::getline(ss_b, item, '|')) parts_b.push_back(item);
        
        if (parts_a.size() >= 4 && parts_b.size() >= 4) {
            // compare by contig_id first
            if (parts_a[1] != parts_b[1]) {
                return parts_a[1] < parts_b[1];
            }
            // then by out_clip
            return std::stoi(parts_a[3]) < std::stoi(parts_b[3]);
        }
        return a < b; // fallback
    });
    
    // assign event IDs
    key_to_event_id.clear();
    for (size_t i = 0; i < sorted_keys.size(); ++i) {
        key_to_event_id[sorted_keys[i]] = "ev" + std::to_string(i + 1);
    }
    
    cout << "assigned " << key_to_event_id.size() << " event IDs" << endl;
}

void RearrangeManager::process_libraries_with_event_ids()
{
    cout << "processing libraries with assigned event IDs..." << endl;
    
    for (const std::string& lib_id : library_ids) {
        cout << "completing processing for library: " << lib_id << endl;
        
        // assign event IDs to this library
        rearrangers[lib_id]->assign_event_ids(key_to_event_id);
        
        // complete aggregation and coverage calculation
        rearrangers[lib_id]->aggregate_events();
        rearrangers[lib_id]->get_coverage();
        
        // write per-library files
        rearrangers[lib_id]->write_to_csv("output/" + lib_id);
        
        cout << "library " << lib_id << " processing completed" << endl;
    }
}

void RearrangeManager::aggregate_across_libraries()
{
    cout << "aggregating events across libraries..." << endl;
    
    // collect data from all libraries
    for (const std::string& lib_id : library_ids) {
        const auto& lib_events = rearrangers[lib_id]->aggregated_events;
        
        for (const AggregatedEvent& event : lib_events) {
            std::string key = event.create_key();
            
            if (rearrangement_data.find(key) != rearrangement_data.end()) {
                // existing event - update totals
                rearrangement_data[key].total_support += event.read_count;
                rearrangement_data[key].total_coverage += event.read_coverage;
                rearrangement_data[key].lib_support[lib_id] = event.read_count;
                rearrangement_data[key].lib_coverage[lib_id] = event.read_coverage;
                rearrangement_data[key].library_count++;
            } else {
                // new event
                SetRearrangementData data;
                data.type = event.type;
                data.contig_id = event.contig_id;
                data.strand = event.strand;
                data.out_clip = event.out_clip;
                data.in_clip = event.in_clip;
                data.element_contig = event.element_contig;
                data.element_strand = event.element_strand;
                data.element_start = event.element_start;
                data.element_end = event.element_end;
                data.total_support = event.read_count;
                data.total_coverage = event.read_coverage;
                data.library_count = 1;
                data.lib_support[lib_id] = event.read_count;
                data.lib_coverage[lib_id] = event.read_coverage;
                
                rearrangement_data[key] = data;
            }
        }
    }
    
    // convert to output rows
    rearrangement_rows.clear();
    for (const auto& pair : rearrangement_data) {
        const SetRearrangementData& data = pair.second;
        
        RearrangementOutputRow row;
        row.event_id = key_to_event_id[pair.first];
        row.type = data.type;
        row.contig_id = data.contig_id;
        row.strand = data.strand;
        row.out_clip = data.out_clip;
        row.in_clip = data.in_clip;
        row.element_contig = data.element_contig;
        row.element_strand = data.element_strand;
        row.element_start = data.element_start;
        row.element_end = data.element_end;
        row.library_count = data.library_count;
        row.total_support = data.total_support;
        row.total_coverage = data.total_coverage;
        row.frequency = data.get_frequency();
        
        rearrangement_rows.push_back(row);
    }
    
    // sort by total support (descending)
    std::sort(rearrangement_rows.begin(), rearrangement_rows.end(),
              [](const RearrangementOutputRow& a, const RearrangementOutputRow& b) {
                  return a.total_support > b.total_support;
              });
    
    cout << "found " << rearrangement_rows.size() << " unique events across all libraries" << endl;
}

void RearrangeManager::write_to_csv(const std::string& ofn_prefix)
{
    cout << "writing cross-library rearrangement results..." << endl;
    
    write_events_file(ofn_prefix);
    write_support_matrix(ofn_prefix);
    write_coverage_matrix(ofn_prefix);
}

void RearrangeManager::write_events_file(const std::string& ofn_prefix)
{
    string filename = ofn_prefix + "_events.tsv";
    cout << "writing events to " << filename << endl;
    
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "error: could not open file " << filename << endl;
        return;
    }
    
    file << "event_id\ttype\tcontig_id\tstrand\tout_clip\tin_clip\telement_contig\telement_strand\telement_start\telement_end\tlibrary_count\ttotal_support\ttotal_coverage\tfrequency\n";
    
    for (const RearrangementOutputRow& row : rearrangement_rows) {
        file << row.event_id << "\t"
             << row.type << "\t"
             << row.contig_id << "\t"
             << row.strand << "\t"
             << row.out_clip << "\t"
             << row.in_clip << "\t"
             << row.element_contig << "\t"
             << row.element_strand << "\t"
             << row.element_start << "\t"
             << row.element_end << "\t"
             << row.library_count << "\t"
             << row.total_support << "\t"
             << row.total_coverage << "\t"
             << row.frequency << "\n";
    }
    
    file.close();
    cout << "wrote " << rearrangement_rows.size() << " events to " << filename << endl;
}

void RearrangeManager::write_support_matrix(const std::string& ofn_prefix)
{
    string filename = ofn_prefix + "_support.tsv";
    cout << "writing support matrix to " << filename << endl;
    
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "error: could not open file " << filename << endl;
        return;
    }
    
    // write header with event details and library columns
    file << "event_id\ttype\tcontig\tclip_out\tclip_in";
    for (const string& lib_id : library_ids) {
        file << "\t" << lib_id;
    }
    file << "\n";
    
    // write data rows
    for (const RearrangementOutputRow& row : rearrangement_rows) {
        file << row.event_id << "\t"
             << row.type << "\t"
             << row.contig_id << "\t"
             << row.out_clip << "\t"
             << row.in_clip;
        
        const SetRearrangementData& data = rearrangement_data[row.create_key()];
        for (const string& lib_id : library_ids) {
            auto it = data.lib_support.find(lib_id);
            file << "\t" << (it != data.lib_support.end() ? it->second : 0);
        }
        file << "\n";
    }
    
    file.close();
    cout << "wrote support matrix to " << filename << endl;
}

void RearrangeManager::write_coverage_matrix(const std::string& ofn_prefix)
{
    string filename = ofn_prefix + "_coverage.tsv";
    cout << "writing coverage matrix to " << filename << endl;
    
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "error: could not open file " << filename << endl;
        return;
    }
    
    // write header with event details and library columns
    file << "event_id\ttype\tcontig\tclip_out\tclip_in";
    for (const string& lib_id : library_ids) {
        file << "\t" << lib_id;
    }
    file << "\n";
    
    // write data rows
    for (const RearrangementOutputRow& row : rearrangement_rows) {
        file << row.event_id << "\t"
             << row.type << "\t"
             << row.contig_id << "\t"
             << row.out_clip << "\t"
             << row.in_clip;
        
        const SetRearrangementData& data = rearrangement_data[row.create_key()];
        for (const string& lib_id : library_ids) {
            auto it = data.lib_coverage.find(lib_id);
            file << "\t" << (it != data.lib_coverage.end() ? it->second : 0);
        }
        file << "\n";
    }
    
    file.close();
    cout << "wrote coverage matrix to " << filename << endl;
}
