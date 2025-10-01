#include "Sequences.h"
#include "RearrangeReadEvent.h"
#include "utils.h"
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <fstream>

using namespace std;

ResolveSeams string_to_resolve_seams(const string& str) {
    string lower_str = str;
    transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);
    
    if (lower_str == "no") {
        return ResolveSeams::NO;
    } else if (lower_str == "reads_only") {
        return ResolveSeams::READS_ONLY;
    } else if (lower_str == "reference_only") {
        return ResolveSeams::REFERENCE_ONLY;
    } else if (lower_str == "complete") {
        return ResolveSeams::COMPLETE;
    } else {
        cerr << "error: invalid resolve_seams value: " << str << endl;
        cerr << "valid values: no, reads_only, reference_only, complete" << endl;
        exit(1);
    }
}

string resolve_seams_to_string(ResolveSeams mode) {
    switch (mode) {
        case ResolveSeams::NO:
            return "no";
        case ResolveSeams::READS_ONLY:
            return "reads_only";
        case ResolveSeams::REFERENCE_ONLY:
            return "reference_only";
        case ResolveSeams::COMPLETE:
            return "complete";
        default:
            return "unknown";
    }
}

// AssemblySequences implementation
AssemblySequences::AssemblySequences() : loaded(false) {}

void AssemblySequences::load_from_file(const string& filename) {
    cout << "loading assembly sequences from " << filename << endl;
    
    sequences.clear();
    unordered_set<string> empty_set; // load all contigs
    read_fasta(filename, empty_set, sequences);
    
    loaded = true;
    cout << "loaded " << sequences.size() << " assembly sequences" << endl;
}

void AssemblySequences::load_from_map(const unordered_map<string, string>& sequence_map) {
    cout << "loading assembly sequences from map" << endl;
    
    sequences = sequence_map;
    loaded = true;
    
    cout << "loaded " << sequences.size() << " assembly sequences" << endl;
}

string AssemblySequences::get_sequence(const string& contig_id) const {
    auto it = sequences.find(contig_id);
    if (it != sequences.end()) {
        return it->second;
    }
    return ""; // return empty string if not found
}

bool AssemblySequences::has_contig(const string& contig_id) const {
    return sequences.find(contig_id) != sequences.end();
}

void AssemblySequences::clear() {
    sequences.clear();
    loaded = false;
}

unordered_set<string> AssemblySequences::get_contig_ids() const {
    unordered_set<string> ids;
    for (const auto& pair : sequences) {
        ids.insert(pair.first);
    }
    return ids;
}

// ReadSequences implementation
ReadSequences::ReadSequences() {}

void ReadSequences::register_lib_file(const string& lib_id, const string& filename) {
    cout << "registering library " << lib_id << " with file " << filename << endl;
    lib_id_to_filename[lib_id] = filename;
}

void ReadSequences::get_sequences(const string& lib_id, const vector<string>& read_ids,
                                 unordered_map<string, string>& result) const {
    cout << "loading " << read_ids.size() << " sequences for library " << lib_id << "..." << endl;
    
    // clear result map
    result.clear();
    
    // check if library is registered
    auto it = lib_id_to_filename.find(lib_id);
    if (it == lib_id_to_filename.end()) {
        cout << "warning: library " << lib_id << " not registered, no sequences loaded" << endl;
        return;
    }
    
    const string& filename = it->second;
    
    // convert vector to set for faster lookup
    unordered_set<string> reads_needed(read_ids.begin(), read_ids.end());
    
    // determine file type and read accordingly
    FileType file_type = get_file_type(filename);
    switch (file_type) {
        case FileType::FASTA:
            read_fasta(filename, reads_needed, result);
            break;
        case FileType::FASTQ:
            read_fastq(filename, reads_needed, result);
            break;
        default:
            cerr << "error: unsupported file type for " << filename << endl;
            exit(1);
    }
    
    cout << "loaded " << result.size() << " sequences" << endl;
    
    // assert all requested reads were found
    for (const string& read_id : read_ids) {
        massert(result.find(read_id) != result.end(),
               "read_id %s not found in library %s file %s", 
               read_id.c_str(), lib_id.c_str(), filename.c_str());
    }
}


void ReadSequences::clear() {
    lib_id_to_filename.clear();
}

bool ReadSequences::has_lib(const string& lib_id) const {
    return lib_id_to_filename.find(lib_id) != lib_id_to_filename.end();
}

unordered_set<string> ReadSequences::get_lib_ids() const {
    unordered_set<string> ids;
    for (const auto& pair : lib_id_to_filename) {
        ids.insert(pair.first);
    }
    return ids;
}

void ReadSequences::get_sequences_with_cache(const string& lib_id, const vector<string>& read_ids,
                                             unordered_map<string, string>& sequences, const string& output_prefix) const
{
    if (!output_prefix.empty()) {
        string cache_filename = output_prefix + "_reads_" + lib_id + ".fasta";
        
        // check if cache file exists before trying to load
        ifstream cache_check(cache_filename);
        if (cache_check.good()) {
            cache_check.close();
            if (load_sequences_from_cache(cache_filename, read_ids, sequences)) {
                cout << "loaded " << sequences.size() << " sequences from cache: " << cache_filename << endl;
                return;
            }
            // if load failed, fall through to reload from original
        }
        
        cout << "cache miss, loading from original files..." << endl;
        get_sequences(lib_id, read_ids, sequences);
        save_sequences_to_cache(cache_filename, sequences);
    } else {
        // no output prefix provided, skip caching
        get_sequences(lib_id, read_ids, sequences);
    }
}

bool ReadSequences::load_sequences_from_cache(const string& cache_filename, const vector<string>& read_ids, 
                                              unordered_map<string, string>& sequences) const
{
    ifstream cache_file(cache_filename);
    if (!cache_file.is_open()) {
        return false; // should not happen since we checked existence, but safety check
    }
    
    sequences.clear();
    string line;
    string current_read_id;
    string current_sequence;
    
    try {
        while (getline(cache_file, line)) {
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                // save previous sequence if we have one
                if (!current_read_id.empty() && !current_sequence.empty()) {
                    sequences[current_read_id] = current_sequence;
                }
                
                // start new sequence
                current_read_id = line.substr(1); // remove '>'
                current_sequence.clear();
            } else {
                // accumulate sequence
                current_sequence += line;
            }
        }
        
        // save last sequence
        if (!current_read_id.empty() && !current_sequence.empty()) {
            sequences[current_read_id] = current_sequence;
        }
        
        cache_file.close();
        
        // verify we have all needed sequences
        for (const string& read_id : read_ids) {
            if (sequences.find(read_id) == sequences.end()) {
                cout << "cache incomplete, missing read: " << read_id << endl;
                sequences.clear();
                return false;
            }
        }
        
        return true;
        
    } catch (const exception& e) {
        cout << "cache file corrupted: " << e.what() << ", deleting and falling back to original" << endl;
        cache_file.close();
        remove(cache_filename.c_str()); // delete corrupted cache
        sequences.clear();
        return false;
    }
}

void ReadSequences::save_sequences_to_cache(const string& cache_filename, const unordered_map<string, string>& sequences) const
{
    if (sequences.empty()) {
        return; // nothing to save
    }
    
    ofstream cache_file(cache_filename);
    if (!cache_file.is_open()) {
        cout << "warning: could not create cache file: " << cache_filename << endl;
        return;
    }
    
    try {
        for (const auto& seq_pair : sequences) {
            const string& read_id = seq_pair.first;
            const string& sequence = seq_pair.second;
            
            cache_file << ">" << read_id << "\n";
            cache_file << sequence << "\n";
        }
        
        cache_file.close();
        cout << "saved " << sequences.size() << " sequences to cache: " << cache_filename << endl;
        
    } catch (const exception& e) {
        cout << "warning: failed to write cache file: " << e.what() << endl;
        cache_file.close();
        remove(cache_filename.c_str()); // delete incomplete cache
    }
}