# cov_intervals

Computes per-segment coverage: for each interval in a segment table, reports the fraction of bases covered by at least one passing alignment (`x_coverage`). Useful for assessing how well a source genome maps to a reference assembly at the bin level.

## Syntax

```bash
alntools cov_intervals -ifn_aln <input.aln> -ifn_segments <segments.txt> -ofn <output.txt> [options]
```

## Parameters

**Mandatory Arguments:**
* `-ifn_aln <fn>`: Input ALN file.
* `-ifn_segments <fn>`: Tab-delimited segment table with at minimum `contig`, `start` (1-based), `end` (1-based inclusive) columns.
* `-ofn <fn>`: Output file path.

**Optional Arguments (alignment filtering):**
* `-clip_mode <str>`: Clip mode for filtering (default: `complete`). Options: `all`, `complete`, `allow_one_side_clip`, `only_one_side_clipped`, `only_two_side_clipped`, `only_clipped`, `local_align`.
* `-clip_margin <int>`: Clip margin in bases (default: `10`).
* `-min_mutations_percent <float>`: Minimum mutation percentage (default: `0.0`).
* `-max_mutations_percent <float>`: Maximum mutation percentage (default: `0.1`).
* `-min_alignment_length <int>`: Minimum alignment length in bases (default: `1000`; use `0` for no limit).
* `-max_alignment_length <int>`: Maximum alignment length in bases (default: `0` = no limit).
* `-min_indel_length <int>`: Minimum indel length counted as a mutation (default: `3`).

## Example

```bash
alntools cov_intervals \
    -ifn_aln output/align.aln \
    -ifn_segments output/focal_segments.txt \
    -ofn output/seg_cov.txt \
    -clip_mode complete \
    -clip_margin 10 \
    -min_mutations_percent 0.0 \
    -max_mutations_percent 0.1 \
    -min_alignment_length 1000
```

## Input File Formats

- **ALN file**: Binary alignment file created by the `construct` command.
- **Segment table**: Tab-delimited file with header. Required columns: `contig`, `start`, `end`. All other columns are passed through unchanged.

```
segment  contig  start  end    length  bin
s1       ctg0    1      14091  14091   b2693
s2       ctg1    1      4193   4193    b5219
```

## Output File Format

The input segment table with one additional column appended:

- `x_coverage`: Fraction of the segment's bases covered by at least one passing alignment (0.0–1.0).

```
segment  contig  start  end    length  bin    x_coverage
s1       ctg0    1      14091  14091   b2693  0.847
s2       ctg1    1      4193   4193    b5219  0.0
```

## Notes

- Only alignments that pass all filtering criteria contribute to coverage.
- Coverage is computed as a boolean bitmap per segment: each base is either covered (1) or not (0), regardless of alignment depth. The `x_coverage` value is the fraction of covered bases.
- Segments whose contig is absent from the ALN file receive `x_coverage = 0`.
- For bin-level summary, aggregate by the `bin` column using a weighted mean: `sum(x_coverage * length) / sum(length)`.
