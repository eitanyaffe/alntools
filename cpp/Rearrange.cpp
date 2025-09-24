#include "Rearrange.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>

using namespace std;

Rearrange::Rearrange(const map<string, AlignmentStore>& stores,
                     const vector<Interval>& intervals,
                     RearrangeVerify* verifier,
                     int max_margin,
                     int min_element_length,
                     int min_anchor_length,
                     double max_anchor_mutations_percent,
                     double max_element_mutation_percent,
                     bool write_per_library_files)
    : stores(stores)
    , intervals(intervals)
    , max_margin(max_margin)
    , min_element_length(min_element_length)
    , min_anchor_length(min_anchor_length)
    , max_anchor_mutations_percent(max_anchor_mutations_percent)
    , max_element_mutation_percent(max_element_mutation_percent)
    , verifier(verifier)
    , write_per_library_files(write_per_library_files)
    , event_grouper(max_margin)
{
    // extract ordered library IDs
    for (const auto& entry : stores) {
        library_ids.push_back(entry.first);
    }
}

void Rearrange::execute()
{
    cout << "executing rearrangement detection across " << library_ids.size() << " libraries..." << endl;
    
    // step 1: detect events per library
    detect_events_per_library();
    
    // step 2: group events and select representatives
    group_and_select_representatives();
    
    // step 3: calculate coverage matrices
    calculate_coverage_matrices();
    
    // step 4: prepare output tables
    prepare_output_tables();
    
    cout << "rearrangement detection completed" << endl;
}

void Rearrange::detect_events_per_library()
{
    cout << "detecting events per library..." << endl;
    
    for (const string& lib_id : library_ids) {
        cout << "--------------------------------" << endl;
        cout << "processing library: " << lib_id << endl;
        
        // create worker on stack
        RearrangeReadEvent detector(stores.at(lib_id), verifier, max_margin,
                                   min_element_length, min_anchor_length,
                                   max_anchor_mutations_percent, max_element_mutation_percent);
        
        // apply interval filtering
        detector.set_intervals(intervals);
        
        // prepare output containers
        vector<ReadEvent>& events = lib_events[lib_id];
        map<string, size_t> rejections;
        
        // detect events in this library
        detector.detect_events(events, rejections);
                
        // print rejection summary for this library
        size_t total_rejected = detector.get_total_rejected_count(rejections);
        
        cout << "library " << lib_id << ": found " << lib_events[lib_id].size() << " events, rejected " 
             << total_rejected << " candidates" << endl;
        
        if (total_rejected > 0) {
            cout << "  rejection breakdown:" << endl;
            for (const auto& pair : rejections) {
                cout << "    " << pair.first << ": " << pair.second << endl;
            }
        }
    }
}

void Rearrange::group_and_select_representatives()
{
    cout << "grouping events and selecting representatives..." << endl;
    
    // perform grouping - pass containers by reference, worker will clear them
    event_grouper.group_events(lib_events, rep_events, event_groups, support_matrix);
    
    cout << "selected " << rep_events.size() << " representative events" << endl;
}

void Rearrange::calculate_coverage_matrices()
{
    cout << "calculating coverage matrices..." << endl;
    
    // calculate coverage per library
    for (const string& lib_id : library_ids) {
        cout << "calculating coverage for library: " << lib_id << endl;
        
        // create worker on stack and calculate coverage directly into matrix
        RearrangeCoverage calculator(stores.at(lib_id));
        calculator.calculate_coverage(rep_events, coverage_matrix[lib_id]);
    }
    
    cout << "coverage calculation completed" << endl;
}

void Rearrange::prepare_output_tables()
{
    cout << "preparing output tables (lib_events size: " << lib_events.size() << ")..." << endl;
    
    // prepare indexed events list for mapping
    vector<tuple<string, size_t, const ReadEvent*>> all_indexed_events;
    for (const auto& lib_pair : lib_events) {
        const string& lib_id = lib_pair.first;
        const vector<ReadEvent>& events = lib_pair.second;
        for (size_t i = 0; i < events.size(); ++i) {
            all_indexed_events.emplace_back(lib_id, i, &events[i]);
        }
    }
    
    // convert to ReadEventRow format using event_groups
    read_event_rows.clear();
    
    for (size_t group_idx = 0; group_idx < event_groups.size(); ++group_idx) {
        const string& event_id = rep_events[group_idx].event_id;
        
        for (int event_idx : event_groups[group_idx]) {
            if (event_idx < static_cast<int>(all_indexed_events.size())) {
                const auto& indexed_event = all_indexed_events[event_idx];
                const string& lib_id = get<0>(indexed_event);
                const ReadEvent& read_event = *get<2>(indexed_event);
                
                ReadEventRow row;
                row.lib_id = lib_id;
                row.read_id = read_event.read_id;
                row.event_id = event_id;
                row.type = read_event.type;
                row.contig_id = read_event.contig_id;
                row.read_strand = read_event.read_strand;
                row.out_clip = read_event.out_clip;
                row.in_clip = read_event.in_clip;
                row.read_clip_out = read_event.read_clip_out;
                row.read_clip_in = read_event.read_clip_in;
                row.element_contig = read_event.element_contig;
                row.element_strand = read_event.element_strand;
                row.element_start = read_event.element_start;
                row.element_end = read_event.element_end;
                row.left_shim = read_event.left_shim;
                row.right_shim = read_event.right_shim;
                row.middle_shim = read_event.middle_shim;
                
                read_event_rows.push_back(row);
            }
        }
    }
    
    cout << "prepared " << read_event_rows.size() << " read event rows and " 
         << rep_events.size() << " event rows" << endl;
}

void Rearrange::load_read_sequences_from_file(const string& filename)
{
    // TODO: implement when needed
    cout << "loading read sequences from file: " << filename << " (not implemented)" << endl;
}

void Rearrange::load_read_sequences_per_library(const map<string, string>& lib_to_read_file)
{
    // TODO: implement when needed
    cout << "loading per-library read sequences (not implemented)" << endl;
}

void Rearrange::load_read_sequences_from_map(const map<string, unordered_map<string, string>>& lib_to_sequences)
{
    // TODO: implement when needed
    cout << "loading read sequences from map (not implemented)" << endl;
}

map<string, size_t> Rearrange::get_rejection_counts() const
{
    map<string, size_t> total_rejections;
    
    // note: we can't collect from detectors here since they're created on the stack
    // rejection counts are printed during execution
    // this method is mainly for R interface compatibility
    
    return total_rejections;
}

void Rearrange::write_to_csv(const string& ofn_prefix)
{
    cout << "writing results to CSV files with prefix: " << ofn_prefix << endl;
    
    write_read_events_file(ofn_prefix);
    write_events_file(ofn_prefix);
    write_support_matrix_file(ofn_prefix);
    write_coverage_matrix_file(ofn_prefix);
    
    cout << "CSV output completed" << endl;
}

void Rearrange::write_read_events_file(const string& ofn_prefix)
{
    string filename = ofn_prefix + "_read_events.txt";
    ofstream ofs(filename);
    
    if (!ofs.is_open()) {
        cerr << "error: cannot open file " << filename << " for writing" << endl;
        return;
    }
    
    // write header
    ofs << "lib_id\tread_id\tevent_id\ttype\tcontig_id\tread_strand\tout_clip\tin_clip"
        << "\tread_clip_out\tread_clip_in\telement_contig\telement_strand"
        << "\telement_start\telement_end\tleft_shim\tright_shim\tmiddle_shim" << endl;
    
    // write data
    for (const ReadEventRow& row : read_event_rows) {
        ofs << row.lib_id << "\t" << row.read_id << "\t" << row.event_id << "\t"
            << row.type << "\t" << row.contig_id << "\t" << row.read_strand << "\t"
            << row.out_clip << "\t" << row.in_clip << "\t" << row.read_clip_out << "\t"
            << row.read_clip_in << "\t" << row.element_contig << "\t" << row.element_strand << "\t"
            << row.element_start << "\t" << row.element_end << "\t" << row.left_shim << "\t"
            << row.right_shim << "\t" << row.middle_shim << endl;
    }
    
    ofs.close();
    cout << "wrote " << read_event_rows.size() << " read events to " << filename << endl;
}

void Rearrange::write_events_file(const string& ofn_prefix)
{
    string filename = ofn_prefix + "_events.txt";
    ofstream ofs(filename);
    
    if (!ofs.is_open()) {
        cerr << "error: cannot open file " << filename << " for writing" << endl;
        return;
    }
    
    // write header
    ofs << "event_id\ttype\tcontig_id\tout_clip\tin_clip\telement_contig"
        << "\telement_strand\telement_start\telement_end\tleft_shim\tright_shim\tmiddle_shim" << endl;
    
    // write data
    for (const Event& event : rep_events) {
        ofs << event.event_id << "\t" << event.type << "\t" << event.contig_id << "\t"
            << event.out_clip << "\t" << event.in_clip << "\t" << event.element_contig << "\t"
            << event.element_strand << "\t" << event.element_start << "\t" << event.element_end << "\t"
            << event.left_shim << "\t" << event.right_shim << "\t" << event.middle_shim << endl;
    }
    
    ofs.close();
    cout << "wrote " << rep_events.size() << " events to " << filename << endl;
}

void Rearrange::write_support_matrix_file(const string& ofn_prefix)
{
    string filename = ofn_prefix + "_support_matrix.txt";
    ofstream ofs(filename);
    
    if (!ofs.is_open()) {
        cerr << "error: cannot open file " << filename << " for writing" << endl;
        return;
    }
    
    // write header
    ofs << "event_id";
    for (const string& lib_id : library_ids) {
        ofs << "\t" << lib_id;
    }
    ofs << endl;
    
    // write data
    for (const Event& event : rep_events) {
        ofs << event.event_id;
        for (const string& lib_id : library_ids) {
            int count = 0;
            auto lib_it = support_matrix.find(lib_id);
            if (lib_it != support_matrix.end()) {
                auto event_it = lib_it->second.find(event.event_id);
                if (event_it != lib_it->second.end()) {
                    count = event_it->second;
                }
            }
            ofs << "\t" << count;
        }
        ofs << endl;
    }
    
    ofs.close();
    cout << "wrote support matrix to " << filename << endl;
}

void Rearrange::write_coverage_matrix_file(const string& ofn_prefix)
{
    string filename = ofn_prefix + "_coverage_matrix.txt";
    ofstream ofs(filename);
    
    if (!ofs.is_open()) {
        cerr << "error: cannot open file " << filename << " for writing" << endl;
        return;
    }
    
    // write header
    ofs << "event_id";
    for (const string& lib_id : library_ids) {
        ofs << "\t" << lib_id;
    }
    ofs << endl;
    
    // write data
    for (const Event& event : rep_events) {
        ofs << event.event_id;
        for (const string& lib_id : library_ids) {
            int count = 0;
            auto lib_it = coverage_matrix.find(lib_id);
            if (lib_it != coverage_matrix.end()) {
                auto event_it = lib_it->second.find(event.event_id);
                if (event_it != lib_it->second.end()) {
                    count = event_it->second;
                }
            }
            ofs << "\t" << count;
        }
        ofs << endl;
    }
    
    ofs.close();
    cout << "wrote coverage matrix to " << filename << endl;
}
