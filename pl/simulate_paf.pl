#!/usr/bin/env perl
use strict;
use warnings;

# parameters for simulation
my $contig_name = "contig1";
my $contig_length = 10000;
my $min_coord = 1;
my $max_coord = 1000;
my $alignment_length = 100;
my $max_mutations = 3;
my $num_alignments = 100;

# output file name
if (@ARGV == 0) {
    print "usage: $0 <output_filename>\n";
    exit 1;
}
my $output_filename = $ARGV[0];

# open output file
open(my $fh, '>', $output_filename) or die "error: could not open output file $output_filename: $!";

# generate alignments
for my $i (1..$num_alignments) {
    my $alignment = generate_alignment(
        $i, $contig_name, $contig_length, $min_coord, $max_coord, 
        $alignment_length, $max_mutations
    );
    print $fh "$alignment\n";
}

close($fh);
print "generated $num_alignments alignments to $output_filename\n";

# function to generate a random alignment with specified constraints
sub generate_alignment {
    my ($alignment_index, $contig_name, $contig_length, $min_coord, $max_coord, 
        $alignment_length, $max_mutations) = @_;
    
    # random start position between min_coord and max_coord-alignment_length
    my $start_pos = int(rand($max_coord - $min_coord - $alignment_length + 1)) + $min_coord;
    my $end_pos = $start_pos + $alignment_length;
    
    # random number of mutations (0 to max_mutations)
    my $num_mutations = int(rand($max_mutations + 1));
    
    # create fields for PAF format
    my $query_name = "read_$alignment_index";
    my $query_length = $alignment_length;
    my $query_start = 0;
    my $query_end = $alignment_length;
    my $strand = '+';
    my $target_length = $contig_length;
    my $matching_bases = $alignment_length - $num_mutations;
    my $alignment_block_length = $alignment_length;
    my $mapping_quality = 60;
    
    # generate CS tag and corresponding CIGAR string
    my ($cs_tag, $cg_tag) = generate_cs_tag($alignment_length, $num_mutations);
    
    # format PAF line
    my $paf_line = join("\t", 
        $query_name,
        $query_length,
        $query_start,
        $query_end,
        $strand,
        $contig_name,
        $target_length,
        $start_pos,
        $end_pos,
        $matching_bases,
        $alignment_block_length,
        $mapping_quality,
        "NM:i:$num_mutations",
        "cg:Z:$cg_tag",
        "cs:Z:$cs_tag"
    );
    
    return $paf_line;
}

# generate CS tag for simulating mutations
sub generate_cs_tag {
    my ($alignment_length, $num_mutations) = @_;
    
    # nucleotide bases
    my @bases = ('a', 'c', 'g', 't');
    
    # if no mutations, return simple match
    if ($num_mutations == 0) {
        return (":$alignment_length", "${alignment_length}M");
    }
    
    # determine mutation positions
    my @positions = ();
    for (my $i = 0; $i < $num_mutations; $i++) {
        my $pos;
        do {
            $pos = int(rand($alignment_length));
        } while (grep { $_ == $pos } @positions);
        push @positions, $pos;
    }
    @positions = sort { $a <=> $b } @positions;
    
    # randomly select mutation types (mismatch, insertion, deletion)
    my @mutation_types = ();
    for (my $i = 0; $i < $num_mutations; $i++) {
        push @mutation_types, int(rand(3)); # 0=mismatch, 1=insertion, 2=deletion
    }
    
    # build CS string
    my $cs_string = "";
    my $cigar_string = "";
    my $last_pos = 0;
    
    for (my $i = 0; $i < $num_mutations; $i++) {
        my $pos = $positions[$i];
        my $type = $mutation_types[$i];
        
        # add match segment before mutation
        if ($pos > $last_pos) {
            $cs_string .= ":" . ($pos - $last_pos);
            $cigar_string .= ($pos - $last_pos) . "M";
        }
        
        # add mutation
        my $ref_base = $bases[int(rand(4))];
        my $alt_base = $bases[int(rand(4))];
        while ($alt_base eq $ref_base) {
            $alt_base = $bases[int(rand(4))];
        }
        
        if ($type == 0) {  # mismatch
            $cs_string .= "*$ref_base$alt_base";
            $cigar_string .= "1M";
            $last_pos = $pos + 1;
        } elsif ($type == 1) {  # insertion
            $cs_string .= "+$alt_base";
            $cigar_string .= "1I";
            $last_pos = $pos;
        } elsif ($type == 2) {  # deletion
            $cs_string .= "-$ref_base";
            $cigar_string .= "1D";
            $last_pos = $pos + 1;
        }
    }
    
    # add final match segment
    if ($last_pos < $alignment_length) {
        $cs_string .= ":" . ($alignment_length - $last_pos);
        $cigar_string .= ($alignment_length - $last_pos) . "M";
    }
    
    return ($cs_string, $cigar_string);
} 