# extract

Saves ALN file to tab-delimited tables for inspection and analysis.

## Syntax

```bash
alntools extract -ifn <input.aln> -ofn_prefix <output_prefix>
```

## Parameters

**Mandatory Arguments:**
* `-ifn <fn>`: Input ALN file.
* `-ofn_prefix <fn>`: Output prefix for result tables.

## Example

```bash
alntools extract -ifn output/test.aln -ofn_prefix output/extracted
```

## Input File Formats

- **ALN file**: Binary alignment file created by the `construct` command

## Output File Formats

The command exports the binary ALN data to human-readable tab-delimited tables. The exact output format matches the internal structure of the ALN file and includes:

- Alignment information (read IDs, contig IDs, coordinates, mutations)
- Read metadata
- Contig information
- Mutation details

See [File Formats](../file-formats.md) for detailed specifications of the output table formats.

## Notes

- This command is useful for debugging and manual inspection of alignment data
- The exported tables can be imported into R, Python, or other analysis tools
- Use this when you need to examine the detailed structure of your ALN file
- The output tables are much larger than the binary ALN format, so use with caution for large datasets
