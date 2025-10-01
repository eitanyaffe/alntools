# info

Provides basic statistics about an ALN file.

## Syntax

```bash
alntools info -ifn <input.aln>
```

## Parameters

**Mandatory Arguments:**
* `-ifn <fn>`: Input ALN file.

## Example

```bash
alntools info -ifn output/test.aln
```

## Input File Formats

- **ALN file**: Binary alignment file created by the `construct` command

## Output Information

The command outputs the following statistics to standard output:

- Total alignments
- Total reads
- Average alignment length
- Total mutations
- Average mutations per alignment

## Example Output

```
Reading alignment file: output/test.aln
Reads: 1000
Alignments: 1250
Average alignment length: 2847.3
Total mutations: 15432
Average mutations per alignment: 12.3
```

## Notes

- This is a quick way to get basic statistics about your alignment data
- Use this command to verify that your ALN file was constructed correctly
- The statistics help assess the quality and characteristics of your alignment dataset
