# alntools File Formats

This document details the file formats used by `alntools` for both input and output.

Positions are 1-based in input and output files, and represented internalls using a 0-based coordinate system.

## Input Format (PAF)

`alntools` processes alignment data in the Pairwise mApping Format (PAF). The following columns are required:

| Column Index | Description                                      | Type  |
|--------------|--------------------------------------------------|-------|
| 1            | Query sequence name                              | string|
| 2            | Query sequence length                            | int   |
| 3            | Query start                                      | int   |
| 4            | Query end                                        | int   |
| 5            | Relative strand ('+' or '-')                     | char  |
| 6            | Target sequence name (contig ID)                 | string|
| 7            | Target sequence length                           | int   |
| 8            | Target start on original strand                  | int   |
| 9            | Target end on original strand                    | int   |

Additionally, the **`cs:Z:` tag** (difference string) **must be present** in one of the optional fields (column 13 onwards). This tag encodes the base differences between the query (read) and the target (contig) and is essential for mutation analysis.

### Example PAF Line

```text
read_298    150 10  145 +   contig_12   5000    1500    1635    120 135 60  cs:Z::10*at:5+gg:20-c:15*cg:70
```

For a complete example of a PAF file, see `examples/align_100.paf` in the repository.

## Intervals File Format

The intervals file is a tab-delimited file specifying regions to query:

| Column  | Description                           | Type  |
|---------|---------------------------------------|-------|
| contig  | Name of the contig                    | string|
| start   | Start position                        | int   |
| end     | End position                          | int   |

Example intervals file:
```
contig  start   end
ctg25860    120448  130902
ctg26175    230065  252386
```

## Query Output Formats

### 1. Full Mode Output

Produces two tab-delimited files:

#### *_alignments.tsv:

| Column         | Description                              | Type    |
|----------------|------------------------------------------|---------|
| alignment_index| Unique index for the alignment           | int     |
| read_id        | ID of the read                           | string  |
| contig_id      | ID of the contig                         | string  |
| read_start     | Start position on read                   | int     |
| read_end       | End position on read                     | int     |
| contig_start   | Start position on contig                 | int     |
| contig_end     | End position on contig                   | int     |
| is_reverse     | Whether alignment is on reverse strand   | boolean |
| cs_tag         | CIGAR string encoding differences        | string  |
| num_mutations  | Number of mutations in the alignment     | int     |
| height         | Vertical position for visualization      | int     |

#### *_mutations.tsv:

| Column         | Description                              | Type    |
|----------------|------------------------------------------|---------|
| alignment_index| Index of parent alignment                | int     |
| read_id        | ID of the read                           | string  |
| contig_id      | ID of the contig                         | string  |
| type           | Mutation type (SUB/INS/DEL)              | string  |
| position       | Position on contig                       | int     |
| desc           | Description of the mutation              | string  |
| height         | Vertical position for visualization      | int     |

### 2. Pileup Mode Output

Produces *_pileup.tsv:

| Column   | Description                                 | Type   |
|----------|---------------------------------------------|--------|
| contig   | Contig ID                                   | string |
| position | Position on contig                          | int    |
| variant  | Variant type (REF or mutation description)  | string |
| count    | Count of occurrences                        | int    |
| coverage | Total read coverage at this position        | int    |
| cumsum   | Cumulative count                            | int    |

### 3. Bin Mode Output

Produces *_bins.tsv:

| Column          | Description                               | Type   |
|-----------------|-------------------------------------------|--------|
| contig          | Contig ID                                 | string |
| start           | Bin start position                        | int    |
| end             | Bin end position                          | int    |
| length          | Bin length                                | int    |
| read_count      | Number of alignments in bin               | int    |
| mutation_count  | Number of mutations in bin                | int    | 
