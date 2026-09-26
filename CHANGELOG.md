# alntools version log

The version is kept in `VERSION` and printed by `alntools version` (also `-v`, `--version`).
The ALN binary format is unchanged since 1.00, so existing `.aln` files do not need to be rebuilt.

## Migrating from 1.00

Rebuild the tool and rerun the commands whose output changed (listed below). In short:

- `segments`: more breakpoints on reverse-strand reads; the `length` column is now one larger.
- `seg_matrix`: adjacency and reach counts change.
- `cov_matrix`: the `var_N` columns change, and from 1.03 the default alignment filter changes too.
- `csegment_coverage`, `get_read_ids`: from 1.03 the default alignment filter changes.
- `query` bin and variants modes: left/right clip counts on reverse-strand alignments are corrected.

To keep the 1.00 alignment filter in `cov_matrix`, `csegment_coverage` and `get_read_ids`, pass
`-clip_mode complete`. The 1.00 `cov_matrix` variance (mean depth divided by segment length) is no
longer available; `-variance_mode poisson` gives variance equal to the mean depth, and the default
`empirical` gives the measured per-base variance. For MetaBAT2 input this matters: the 1.00 values were
far below 1, and MetaBAT2 raises any variance below 1 to 1.

## 1.03

Defaults now match the settings used for long-read binning.

- `cov_matrix`: default `-clip_mode end_unique` (was `complete`) and default `-variance_mode empirical`
  (was `poisson`).
- `csegment_coverage`, `get_read_ids`: default `-clip_mode end_unique` (was `complete`).
- `cov_intervals`: the `local_align` and `end_unique` clip modes now work (the read index was not built,
  so these modes stopped with an error).
- Docs: `cov_matrix` algorithm corrected (coverage is aligned bases divided by segment length, not a
  read count); library tables use header columns `lib_id` and `aln_fn` (some docs said `id fn`); README
  gains an alignment filtering section.

Output changes: calls that did not set `-clip_mode` or `-variance_mode` in the commands above.
Calls that set both explicitly give the same output as 1.02.

## 1.02

New options only; defaults unchanged, so output is the same as 1.01.

- New clip mode `end_unique` in `cov_matrix`, `csegment_coverage` and `get_read_ids`.
  It accepts every `complete` alignment, plus an alignment that leaves one read end unaligned where the
  alignment reaches the contig end (a read overhanging the contig), provided no other alignment of the
  same read overlaps more than half of it on the read. Alignments that leave a read end unaligned inside
  the contig, or both read ends unaligned, are rejected. With `complete`, long-read depth drops toward
  contig ends over about one read length, because reads running past the contig end are discarded;
  `end_unique` recovers that depth without counting reads at secondary placements.
- New `cov_matrix -variance_mode poisson|empirical` (default `poisson` in 1.02). `empirical` reports the
  population variance of the per-base depth along the segment, computed by sweeping the start and end
  positions of the counted alignments. The mean depth is identical in both modes.

## 1.01

Bug fixes. Output changes for the commands listed.

- `segments`: the test for whether a read simply continues along the contig now follows the strand of
  the anchor alignment. Before, reverse-strand reads could be taken to break where they continued, and
  vice versa, so breakpoints were missed or called in error. A read whose next alignment slightly
  overlaps the anchor on the read is now also treated as continuing. Breakpoints at contig coordinate 0
  are no longer skipped. The `length` column of the segment table is now `end - start + 1` (it was one
  short).
- `seg_matrix`: alignments chosen per read are processed in read order (they were in selection order),
  so adjacency and reach counts change.
- `cov_matrix`: the aligned-base total is 64-bit (a 32-bit total could overflow on very deep segments),
  and `var_N` is the mean depth (it was the mean depth divided by segment length, close to 0).
- `query` bin and variants modes: left and right clipping are assigned by strand; on reverse-strand
  alignments the two were swapped. Variants mode counts an alignment once when it spans several query
  intervals (it was counted once per interval).
- Alignment filter: fixed an unsigned underflow in the test for an alignment reaching the read end,
  which affected reads shorter than `-clip_margin`.
- Mapping a contig position to a read position (used by `seg_matrix` and `get_local_deletions`) is
  corrected for positions inside a deletion.

## 1.00

Baseline release with version reporting (`alntools version`).
