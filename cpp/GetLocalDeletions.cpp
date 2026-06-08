#include "GetLocalDeletions.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

GetLocalDeletions::GetLocalDeletions()
    : margin(10)
    , flank(100)
    , min_indel_length(3)
{
}

string GetLocalDeletions::format_flanking_mutations(const Alignment& aln,
                                                     const AlignmentStore& store,
                                                     uint32_t range_start,
                                                     uint32_t range_end) const
{
    string result;
    for (uint32_t mut_idx : aln.mutations) {
        if (should_skip_short_indel(store, aln.contig_index, mut_idx))
            continue;
        const Mutation& mut = store.get_mutation(aln.contig_index, mut_idx);
        if (mut.position < range_start || mut.position > range_end)
            continue;
        if (!result.empty())
            result += ",";
        result += to_string(mut.position) + "[" + mut.to_string() + "]";
    }
    return result.empty() ? "." : result;
}

void GetLocalDeletions::find_deletions(const AlignmentStore& store,
                                        const vector<Interval>& intervals,
                                        ofstream& out,
                                        bool emit_interval_id) const
{
    int count = 0;
    for (const Interval& interval : intervals) {
        uint32_t q_start = interval.start;
        uint32_t q_end   = interval.end;

        auto alignments = store.get_alignments_intersecting_interval(interval);
        for (const auto& aln_ref : alignments) {
            const Alignment& aln = aln_ref.get();

            for (uint32_t mut_idx : aln.mutations) {
                const Mutation& mut = store.get_mutation(aln.contig_index, mut_idx);
                if (mut.type != MutationType::DELETION)
                    continue;

                uint32_t del_start = mut.position;
                uint32_t del_end   = del_start + static_cast<uint32_t>(mut.nts.size());

                // deletion must approximately match both endpoints of the query interval
                int64_t start_diff = static_cast<int64_t>(del_start) - static_cast<int64_t>(q_start);
                int64_t end_diff   = static_cast<int64_t>(del_end)   - static_cast<int64_t>(q_end);
                if (llabs(start_diff) >= margin || llabs(end_diff) >= margin)
                    continue;

                uint32_t read_position = contig_to_read_coord(aln, del_start, store);

                char aln_strand = aln.is_reverse ? '-' : '+';
                // XOR: same sign → element is in '+' orientation in the read, opposite → '-'
                char del_orientation = ((interval.strand == '+') != aln.is_reverse) ? '+' : '-';

                // pre-window: [del_start - flank, del_start - 1]
                uint32_t pre_start = (del_start >= static_cast<uint32_t>(flank))
                                     ? del_start - static_cast<uint32_t>(flank) : 0;
                string pre  = format_flanking_mutations(aln, store, pre_start, del_start - 1);

                // post-window: [del_end, del_end + flank]
                string post = format_flanking_mutations(aln, store, del_end,
                                                        del_end + static_cast<uint32_t>(flank));

                out << store.get_read_id(aln.read_index)    << "\t"
                    << store.get_contig_id(aln.contig_index) << "\t"
                    << del_start                              << "\t"
                    << del_end                                << "\t"
                    << mut.nts.size()                         << "\t"
                    << aln_strand                             << "\t"
                    << del_orientation                        << "\t"
                    << read_position                          << "\t"
                    << pre                                    << "\t"
                    << post;
                if (emit_interval_id)
                    out << "\t" << interval.name;
                out << "\n";
                count++;
            }
        }
    }
    cout << "found " << count << " matching deletions" << endl;
}

void GetLocalDeletions::compute(const string& ifn_aln,
                                 const string& ifn_intervals,
                                 const string& ofn,
                                 int margin_param,
                                 int flank_param,
                                 int min_indel_length_param)
{
    margin           = margin_param;
    flank            = flank_param;
    min_indel_length = min_indel_length_param;

    vector<Interval> intervals;
    cout << "reading intervals: " << ifn_intervals << endl;
    read_intervals(ifn_intervals, intervals);
    cout << "loaded " << intervals.size() << " intervals" << endl;

    cout << "loading alignment file: " << ifn_aln << endl;
    AlignmentStore store;
    store.load(ifn_aln);

    store.count_short_indels(min_indel_length);

    // emit interval_id column when any interval carries a name
    bool emit_interval_id = false;
    for (const auto& iv : intervals) {
        if (!iv.name.empty()) { emit_interval_id = true; break; }
    }

    ofstream out(ofn.c_str());
    massert(out.is_open(), "could not open output file %s", ofn.c_str());

    out << "read_id\tcontig\tdel_start\tdel_end\tdel_len\taln_strand\tdel_orientation\tread_position\tpre_mutations\tpost_mutations";
    if (emit_interval_id)
        out << "\tinterval_id";
    out << "\n";

    cout << "writing output: " << ofn << endl;
    find_deletions(store, intervals, out, emit_interval_id);
    out.close();

    cout << "get_local_deletions complete" << endl;
}
