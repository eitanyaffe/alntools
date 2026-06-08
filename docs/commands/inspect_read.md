# inspect_read

Prints a tab-delimited dump of all alignments and mutations for a single read in an ALN file. The output is structured (three header-prefixed TSV blocks) so it can be parsed or sliced with standard shell tools. For free-text per-read inspection of one or more reads, use `info -read_ids`.

## Syntax

```bash
alntools inspect_read -ifn_aln <input.aln> -read_id <read_id>
```

## Parameters

**Mandatory Arguments:**
* `-ifn_aln <fn>`: Input ALN file.
* `-read_id <str>`: Read ID to inspect. The command exits with an error if the read ID is not present in the ALN file.

## Input File Formats

- **ALN file**: Binary alignment file created by the `construct` command.

## Output

Written to stdout in three sections, each preceded by a header line.

**Section 1 — Read Summary (free text):**
- Read ID
- Length in bp
- Number of alignments
- Total mutations across all alignments

**Section 2 — Alignments (TSV):**
Columns: `alignment_index`, `contig`, `read_start`, `read_end`, `contig_start`, `contig_end`, `is_reverse`, `mutations`

**Section 3 — Mutations (TSV):**
Columns: `alignment_index`, `contig`, `mutation_index`, `type`, `position`, `description`

`type` is one of `SUB`, `INS`, `DEL`. `description` is the full mutation string (e.g. `A:G`, `+CGT`, `-ATTC`).

## Example

```bash
alntools inspect_read -ifn_aln output/test.aln -read_id read_abc123
```

```
=== Read Summary ===
Read ID: read_abc123
Length: 15420 bp
Number of alignments: 2
Total mutations: 3

=== Alignments ===
alignment_index	contig	read_start	read_end	contig_start	contig_end	is_reverse	mutations
4	contig_1	0	14800	10000	24800	false	3
5	contig_3	14900	15420	5000	5520	true	0

=== Mutations ===
alignment_index	contig	mutation_index	type	position	description
4	contig_1	17	SUB	10452	A:G
4	contig_1	42	DEL	18300	-ATT
4	contig_1	58	INS	22100	+C
```
