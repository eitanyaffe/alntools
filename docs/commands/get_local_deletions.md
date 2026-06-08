# get_local_deletions

Finds reads that carry a deletion closely matching a set of query intervals. For each query interval `[q_start, q_end)`, the command scans alignments overlapping the interval and reports reads whose alignment contains a deletion where both endpoints fall within `margin` bases of the query interval endpoints. The output includes the read ID, the deletion's reference coordinates, the corresponding position in the read, and a summary of nearby mutations flanking the deletion site.

This is useful for identifying reads that lack a mobile element (or other large feature) present in the reference assembly: such reads show a deletion spanning the element's coordinates when aligned to the assembly.

## Syntax

```bash
alntools get_local_deletions -ifn_aln <alignments.aln> \
                              -ifn_intervals <intervals.txt> \
                              -ofn <output.txt> \
                              [-margin <int>] \
                              [-flank <int>] \
                              [-min_indel_length <int>]
```

## Parameters

**Mandatory Arguments:**
* `-ifn_aln <fn>`: Input ALN file.
* `-ifn_intervals <fn>`: Query intervals file (`contig`, `start`, `end` columns).
* `-ofn <fn>`: Output tab-delimited file.

**Optional Arguments:**
* `-margin <int>`: Maximum allowed coordinate difference (in bases) between a candidate deletion endpoint and the corresponding query interval endpoint (default: 10). Both `|del_start - q_start|` and `|del_end - q_end|` must be below this threshold for a deletion to match.
* `-flank <int>`: Window size in bases around the deletion for reporting nearby mutations (default: 100). Mutations in `[del_start - flank, del_start)` are reported as `pre_mutations`; mutations in `[del_end, del_end + flank]` are reported as `post_mutations`.
* `-min_indel_length <int>`: Minimum indel length to include in flanking mutation reporting (default: 3). Shorter indels are excluded from `pre_mutations` and `post_mutations`.

## Examples

**Basic usage:**

```bash
alntools get_local_deletions -ifn_aln output/sample.aln \
                              -ifn_intervals focal_segments.txt \
                              -ofn local_deletions.txt
```

**Tighter matching with wider flanking window:**

```bash
alntools get_local_deletions -ifn_aln output/sample.aln \
                              -ifn_intervals focal_segments.txt \
                              -ofn local_deletions.txt \
                              -margin 5 \
                              -flank 200
```

## Input File Formats

**Intervals file (`ifn_intervals`):**

Tab-delimited file with header. Coordinates are 0-based half-open `[start, end)`. Any accepted header may optionally be followed by a `name` column (e.g. `\tname`). When `name` is present, an `interval_id` column is added to the output.

Accepted base headers (with or without trailing `\tname`):

| Header | Description |
|--------|-------------|
| `contig\tstart\tend` | Basic 3-column format; `del_orientation` will always be `+` |
| `contig\tstart\tend\tstrand` | 4-column format; `strand` is used to compute `del_orientation` |
| `contig\tstart\tend\tvstart\tvend` | With view-coordinates (no strand) |
| `contig\tstart\tend\tvstart\tvend\tstrand` | Full format with view-coordinates and strand |

```
contig	start	end	strand	name
ctg5	412800	421050	+	c4
ctg5	530000	538200	-	c7
```

## Output File Format

Tab-delimited file with header. One row per matching deletion (a single read may contribute multiple rows if it has more than one matching deletion):

| Column | Description |
|--------|-------------|
| `read_id` | Read identifier |
| `contig` | Contig the deletion is on |
| `del_start` | Deletion start on the contig (reference coordinate) |
| `del_end` | Deletion end on the contig (exclusive, reference coordinate) |
| `del_len` | Length of the deletion in bases |
| `aln_strand` | Strand of the read's alignment to the contig (`+` or `-`) |
| `del_orientation` | Orientation of the deleted element within the read (`+` or `-`). Computed as `aln_strand` XOR csegment strand from the interval file. When the interval file has no strand column the interval strand defaults to `+`. |
| `read_position` | Read coordinate of the deletion site. The read bases on either side of the deletion are adjacent in the read, so a single coordinate identifies the gap. |
| `pre_mutations` | Non-short-indel mutations in `[del_start - flank, del_start)`, formatted as `pos[desc]` comma-separated, or `.` if none |
| `post_mutations` | Non-short-indel mutations in `[del_end, del_end + flank]`, formatted as `pos[desc]` comma-separated, or `.` if none |
| `interval_id` | Name of the matched interval (only present when the intervals file contains a `name` column) |

**Mutation descriptor format** (used in `pre_mutations` and `post_mutations`):
- Substitution: `pos[read_base:ref_base]` e.g. `412790[A:T]`
- Insertion: `pos[+seq]` e.g. `412795[+GCA]`
- Deletion: `pos[-seq]` e.g. `412810[-acgt]`

**Example output:**

```
read_id    contig  del_start  del_end  del_len  aln_strand  del_orientation  read_position  pre_mutations  post_mutations
read_00145 ctg5    412803     421047   8244     +           +                31200          412795[A:T]    .
read_00289 ctg5    412801     421049   8248     -           -                18450          .              421052[+GCA]
read_00413 ctg5    412802     421048   8246     +           +                55321          .              .
```

## Algorithm

1. Load the intervals file and ALN file.
2. Call `count_short_indels` with `min_indel_length` to mark short indels in all alignments.
3. For each query interval `[q_start, q_end)`:
   a. Retrieve all alignments intersecting the interval.
   b. For each alignment, scan its mutations for `DELETION` type.
   c. A deletion at `[del_start, del_end)` matches if `|del_start - q_start| < margin` and `|del_end - q_end| < margin`.
   d. For each matching deletion, convert `del_start` to a read coordinate via `contig_to_read_coord`.
   e. Determine `aln_strand` from `alignment.is_reverse` and `del_orientation` from `aln_strand` XOR the interval's strand field.
   f. Collect non-short-indel mutations within `flank` bases on each side.
   g. Emit one output row.

## Notes

- `read_position` reports a single read coordinate at the deletion site. Mathematically, the read coordinates of `del_start` and `del_end` are identical because the deletion's own length cancels in `contig_to_read_coord`.
- `del_orientation` indicates which strand of the deleted element is readable in the original read. `+` means the element's sequence appears in forward orientation in the read; `-` means reverse-complement. This is computed as `(interval_strand == '+') XOR aln.is_reverse`.
- If the interval file has no strand column the interval strand defaults to `+`, so `del_orientation` equals `aln_strand` in that case.
- Short indels (length < `min_indel_length`) are excluded from `pre_mutations` and `post_mutations` but are not excluded from the deletion matching itself; a large deletion is matched regardless of adjacent short indels.
- If a query interval has no matching deletions in any read, it produces no output rows (not a silent zero row).
- Multiple query intervals in a single run are each processed independently; a read matching two different intervals produces one row per interval.
