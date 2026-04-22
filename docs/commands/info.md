# info

Provides basic statistics about an ALN file. When read IDs are provided, prints detailed per-read information instead.

## Syntax

```bash
alntools info -ifn <input.aln> [-read_ids <id1,id2,...>]
```

## Parameters

**Mandatory Arguments:**
* `-ifn <fn>`: Input ALN file.

**Optional Arguments:**
* `-read_ids <str>`: Comma-separated list of read IDs. When provided, prints detailed information for each read instead of file-level statistics.

## Examples

**File statistics:**
```bash
alntools info -ifn output/test.aln
```

**Inspect specific reads:**
```bash
alntools info -ifn output/test.aln -read_ids read1,read2,read3
```

## Input File Formats

- **ALN file**: Binary alignment file created by the `construct` command

## Output

**Without `-read_ids`:** File-level summary statistics printed to stdout:
- Total alignments
- Total reads
- Average alignment length
- Total mutations
- Average mutations per alignment

**With `-read_ids`:** One block per read, separated by `============================================================`. Each block contains:
- Read ID and length
- Number of alignments
- Per-alignment details: contig, contig coordinates, read coordinates, strand, mutation count
- Per-mutation details: type, contig position, description

## Example Output

**Statistics mode:**
```
Total alignments: 1250
Total reads: 1000
Average alignment length: 2847.3 bp
Total mutations: 15432
Average mutations per alignment: 12.3
```

**Read inspection mode:**
```
============================================================
read: read_abc123
length: 15420 bp
alignment count: 2

alignment 1:
  contig: contig_1
  contig coords: 10000-25000
  read coords: 0-14800
  strand: +
  mutation count: 3
  mutation 1: SUB pos=10452 A:G
  mutation 2: DEL pos=18300 -ATT
  mutation 3: INS pos=22100 +C

alignment 2:
  contig: contig_3
  contig coords: 5000-5800
  read coords: 14900-15420
  strand: -
  mutation count: 0
============================================================
```
