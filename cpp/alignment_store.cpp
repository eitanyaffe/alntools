#include "alignment_store.h"
#include "utils.h"
#include <fstream>
#include <iostream>

using std::ofstream;
using std::ifstream;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::unordered_map;
using std::ios;

void AlignmentStore::save(const string& filename) 
{
    ofstream file(filename, ios::binary);
    if (!file.is_open()) {
        cerr << "Error opening file for writing: " << filename << endl;
        abort();
    }

    // Save reads
    size_t num_reads = reads_.size();
    file.write(reinterpret_cast<const char*>(&num_reads), sizeof(num_reads));
    for (const auto& read : reads_) {
        // Write the read id length and id
        size_t id_length = read.id.size();
        file.write(reinterpret_cast<const char*>(&id_length), sizeof(id_length));
        file.write(read.id.c_str(), id_length);

        // Write the read length
        file.write(reinterpret_cast<const char*>(&read.length), sizeof(read.length));
    }
    
    // Save contigs
    size_t num_contigs = contigs_.size();
    file.write(reinterpret_cast<const char*>(&num_contigs), sizeof(num_contigs));
    for (const auto& contig : contigs_) {
        // Write the contig id length and id
        size_t id_length = contig.id.size();
        file.write(reinterpret_cast<const char*>(&id_length), sizeof(id_length));
        file.write(contig.id.c_str(), id_length);

        // Write the contig length
        file.write(reinterpret_cast<const char*>(&contig.length), sizeof(contig.length));
    }

    // Save alignments
    size_t num_alignments = alignments_.size();
    file.write(reinterpret_cast<const char*>(&num_alignments), sizeof(num_alignments));
    for (const auto& alignment : alignments_) {
        // Write the alignment data
        file.write(reinterpret_cast<const char*>(&alignment.read_index), sizeof(alignment.read_index));
        file.write(reinterpret_cast<const char*>(&alignment.contig_index), sizeof(alignment.contig_index));
        file.write(reinterpret_cast<const char*>(&alignment.read_start), sizeof(alignment.read_start));
        file.write(reinterpret_cast<const char*>(&alignment.read_end), sizeof(alignment.read_end));

        // Write mutations
        size_t num_mutations = alignment.mutations.size();
        file.write(reinterpret_cast<const char*>(&num_mutations), sizeof(num_mutations));
        for (const auto& mutation : alignment.mutations) {
            file.write(reinterpret_cast<const char*>(&mutation.type), sizeof(mutation.type));
            file.write(reinterpret_cast<const char*>(&mutation.position), sizeof(mutation.position));

            // Write query_bases
            size_t query_bases_length = mutation.query_bases.size();
            file.write(reinterpret_cast<const char*>(&query_bases_length), sizeof(query_bases_length));
            file.write(mutation.query_bases.c_str(), query_bases_length);

            // Write target_bases
            size_t target_bases_length = mutation.target_bases.size();
            file.write(reinterpret_cast<const char*>(&target_bases_length), sizeof(target_bases_length));
            file.write(mutation.target_bases.c_str(), target_bases_length);
        }
    }

    file.close();
}

void AlignmentStore::load(const string& filename) 
{
    ifstream file(filename, ios::binary);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        abort();
    }

    // Clear existing data
    reads_.clear();
    contigs_.clear();
    read_id_to_index.clear();
    contig_id_to_index.clear();

    // reads
    size_t num_reads = 1;
    file.read(reinterpret_cast<char*>(&num_reads), sizeof(num_reads));
    for (size_t i = 0; i < num_reads; ++i) {
        Read read;

        // Read the read id length and id
        size_t id_length;
        file.read(reinterpret_cast<char*>(&id_length), sizeof(id_length));
        read.id.resize(id_length);
        file.read(&read.id[0], id_length);

        // Read the read length
        file.read(reinterpret_cast<char*>(&read.length), sizeof(read.length));

        reads_.push_back(read);
        read_id_to_index[read.id] = i;
    }

    // Load contigs
    size_t num_contigs;
    file.read(reinterpret_cast<char*>(&num_contigs), sizeof(num_contigs));
    for (size_t i = 0; i < num_contigs; ++i) {
        Contig contig;

        // Read the contig id length and id
        size_t id_length;
        file.read(reinterpret_cast<char*>(&id_length), sizeof(id_length));
        contig.id.resize(id_length);
        file.read(&contig.id[0], id_length);

        // Read the contig length
        file.read(reinterpret_cast<char*>(&contig.length), sizeof(contig.length));

        contigs_.push_back(contig);
        contig_id_to_index[contig.id] = i;
    }

    // Load alignments
    size_t num_alignments;
    file.read(reinterpret_cast<char*>(&num_alignments), sizeof(num_alignments));
    for (size_t i = 0; i < num_alignments; ++i) {
        Alignment alignment;

        // Read the alignment data
        file.read(reinterpret_cast<char*>(&alignment.read_index), sizeof(alignment.read_index));
        file.read(reinterpret_cast<char*>(&alignment.contig_index), sizeof(alignment.contig_index));
        file.read(reinterpret_cast<char*>(&alignment.read_start), sizeof(alignment.read_start));
        file.read(reinterpret_cast<char*>(&alignment.read_end), sizeof(alignment.read_end));

        // Read mutations
        size_t num_mutations;
        file.read(reinterpret_cast<char*>(&num_mutations), sizeof(num_mutations));
        for (size_t j = 0; j < num_mutations; ++j) {
            MutationType type;
            uint32_t position;
            string query_bases;
            string target_bases;
            
            file.read(reinterpret_cast<char*>(&type), sizeof(type));
            file.read(reinterpret_cast<char*>(&position), sizeof(position));
            
            // Read query_bases
            size_t query_bases_length;
            file.read(reinterpret_cast<char*>(&query_bases_length), sizeof(query_bases_length));
            query_bases.resize(query_bases_length);
            file.read(&query_bases[0], query_bases_length);
            
            // Read target_bases
            size_t target_bases_length;
            file.read(reinterpret_cast<char*>(&target_bases_length), sizeof(target_bases_length));
            target_bases.resize(target_bases_length);
            file.read(&target_bases[0], target_bases_length);
            
            // Create mutation with the read data
            Mutation mutation(type, position, query_bases, target_bases);
            alignment.mutations.push_back(mutation);
        }

        alignments_.push_back(alignment);
    }

    file.close();
}

void AlignmentStore::export_tab_delimited(const string& prefix) 
{
    // Write main alignments file
    string alignments_file = prefix + "_alignments.txt";
    std::cout << "Writing alignments to: " << alignments_file << std::endl;
    ofstream alignments_out(alignments_file);
    massert(alignments_out.is_open(), "Failed to open file for writing: %s", alignments_file.c_str());

    // Write mutations file
    string mutations_file = prefix + "_mutations.txt";
    std::cout << "Writing mutations to: " << mutations_file << std::endl;
    ofstream mutations_out(mutations_file);
    massert(mutations_out.is_open(), "Failed to open file for writing: %s", mutations_file.c_str());

    // Write headers
    alignments_out << "read_id\tread_start\tread_end\tcontig_id\tcontig_start\tcontig_end\tmutation_count\n";
    mutations_out << "read_id\tcontig_id\tmutation_type\tposition\tquery_bases\ttarget_bases\n";

    // Write alignments and mutations
    for (const auto& alignment : alignments_) {
        // Get read and contig IDs
        const string& read_id = get_read_id(alignment.read_index);
        const string& contig_id = get_contig_id(alignment.contig_index);

        // Write alignment data
        alignments_out << read_id << "\t"
                      << alignment.read_start << "\t"
                      << alignment.read_end << "\t"
                      << contig_id << "\t"
                      << alignment.contig_start << "\t"
                      << alignment.contig_end << "\t"
                      << alignment.mutations.size() << "\n";

        // Write detailed mutation data
        for (const auto& mutation : alignment.mutations) {
            string mutation_type;
            switch (mutation.type) {
                case MutationType::SUBSTITUTION:
                    mutation_type = "SUB";
                    break;
                case MutationType::INSERTION:
                    mutation_type = "INS";
                    break;
                case MutationType::DELETION:
                    mutation_type = "DEL";
                    break;
                default:
                    massert(false, "Unknown mutation type");
                    break;
            }

            mutations_out << read_id << "\t"
                         << contig_id << "\t"
                         << mutation_type << "\t"
                         << mutation.position << "\t"
                         << mutation.query_bases << "\t"
                         << mutation.target_bases << "\n";
        }
    }

    alignments_out.close();
    mutations_out.close();
}

size_t AlignmentStore::get_read_index(const string& read_id)
{
    auto it = read_id_to_index.find(read_id);
    massert(it != read_id_to_index.end(), "read not found: %s", read_id.c_str());
    return it->second;
}

size_t AlignmentStore::get_contig_index(const string& contig_id)
{
    auto it = contig_id_to_index.find(contig_id);
    massert(it != contig_id_to_index.end(), "contig not found: %s", contig_id.c_str());
    return it->second;
}

size_t AlignmentStore::add_or_get_read_index(const string& read_id, uint32_t length)
{
    auto it = read_id_to_index.find(read_id);
    if (it != read_id_to_index.end()) {
        // Read already exists, return its index
        return it->second;
    } else {
        // Add new Read
        size_t new_index = reads_.size();
        reads_.emplace_back(read_id, length);
        read_id_to_index[read_id] = new_index;
        return new_index;
    }
}

size_t AlignmentStore::add_or_get_contig_index(const string& contig_id, uint32_t length)
{
    auto it = contig_id_to_index.find(contig_id);
    if (it != contig_id_to_index.end()) {
        // Contig already exists, return its index
        return it->second;
    } else {
        // Add new Contig
        size_t new_index = contigs_.size();
        contigs_.emplace_back(contig_id, length);
        contig_id_to_index[contig_id] = new_index;
        return new_index;
    }
}

const string& AlignmentStore::get_read_id(size_t read_index) const
{
    massert(read_index < reads_.size(), "read index out of bounds: %zu", read_index);
    return reads_[read_index].id;
}

const string& AlignmentStore::get_contig_id(size_t contig_index) const
{
    massert(contig_index < contigs_.size(), "contig index out of bounds: %zu", contig_index);
    return contigs_[contig_index].id;
} 
