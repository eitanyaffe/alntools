# Test rules for alntools
# This file contains test rules that are included in the main Makefile

# input
TEST_PAF = examples/align_100.paf
TEST_READS = examples/reads_100.fq
TEST_CONTIGS = examples/contigs_100.fa

# output
TEST_OUTPUT_DIR = output

# intervals for query_full
TEST_INTERVALS_SMALL = examples/intervals_small.txt
TEST_INTERVALS_LARGE = examples/intervals_large.txt

# bin size for query_bin
TEST_BIN_SIZE = 1000

# consensus threshold for query_consensus
TEST_CONSENSUS_THRESHOLD = 0.8
TEST_MIN_CONSENSUS_COVERAGE = 3

.PHONY: test test_basic test_full test_query_full test_query_bin \
test_query_pileup test_query_consensus test_query_variants test_coverage test_query_all test_R_all test_R_commands test_R_plot \
test_create_dense_paf test_genes test_rearrange test_rearrange_basic test_rearrange_reference clean-test test-r-load

########################################################################################
# creating aln files
########################################################################################

# construct ALN w/o validation
test_basic: $(TARGET)
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running BASIC TEST, without paf validation"
	rm -rf $(TEST_OUTPUT_DIR) && mkdir -p $(TEST_OUTPUT_DIR)
	$(TARGET) construct \
		-ifn_paf $(TEST_PAF) \
		-ofn $(TEST_OUTPUT_DIR)/test.aln
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	$(TARGET) info \
		-ifn $(TEST_OUTPUT_DIR)/test.aln
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	$(TARGET) extract \
		-ifn $(TEST_OUTPUT_DIR)/test.aln \
		-ofn_prefix $(TEST_OUTPUT_DIR)/test
	@echo "BASIC TEST completed successfully"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

# construct ALN with validation of cs tags (rarely used)
test_full: $(TARGET)
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running FULL TEST, with paf validation"
	$(TARGET) construct \
		-ifn_paf $(TEST_PAF) \
		-ifn_reads $(TEST_READS) \
		-ifn_contigs $(TEST_CONTIGS) \
		-verify T \
		-ofn $(TEST_OUTPUT_DIR)/test_full.aln
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

########################################################################################
# query tests
########################################################################################

test_query_full: $(TARGET)
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running QUERY FULL"
	$(TARGET) query \
		-ifn_aln $(TEST_OUTPUT_DIR)/test.aln \
		-ifn_intervals $(TEST_INTERVALS_LARGE) \
		-ofn_prefix $(TEST_OUTPUT_DIR)/query \
		-mode full
	@echo "QUERY FULL completed successfully"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

test_query_bin: $(TARGET)
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running QUERY BIN"
	$(TARGET) query \
		-ifn_aln $(TEST_OUTPUT_DIR)/test.aln \
		-ifn_intervals $(TEST_INTERVALS_SMALL) \
		-ofn_prefix $(TEST_OUTPUT_DIR)/query \
		-mode bin \
		-binsize $(TEST_BIN_SIZE)
	@echo "QUERY BIN completed successfully"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

test_query_pileup: $(TARGET)
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running QUERY PILEUP"
	$(TARGET) query \
		-ifn_aln $(TEST_OUTPUT_DIR)/test.aln \
		-ifn_intervals $(TEST_INTERVALS_SMALL) \
		-ofn_prefix $(TEST_OUTPUT_DIR)/query \
		-mode pileup \
		-pileup_mode mutated
	@echo "QUERY PILEUP completed successfully"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

test_query_consensus: $(TARGET)
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running QUERY CONSENSUS"
	$(TARGET) query \
		-ifn_aln $(TEST_OUTPUT_DIR)/test.aln \
		-ifn_intervals $(TEST_INTERVALS_SMALL) \
		-ofn_prefix $(TEST_OUTPUT_DIR)/query \
		-mode consensus \
		-consensus_threshold $(TEST_CONSENSUS_THRESHOLD) \
		-min_consensus_coverage $(TEST_MIN_CONSENSUS_COVERAGE) \
		-num_threads 2
	@echo "QUERY CONSENSUS completed successfully"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

test_query_variants: $(TARGET)
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running QUERY VARIANTS"
	# create libraries table for testing
	@echo "lib_id	aln_fn" > $(TEST_OUTPUT_DIR)/libraries.tsv
	@echo "lib1	$(TEST_OUTPUT_DIR)/test.aln" >> $(TEST_OUTPUT_DIR)/libraries.tsv
	@echo "lib2	$(TEST_OUTPUT_DIR)/test.aln" >> $(TEST_OUTPUT_DIR)/libraries.tsv
	$(TARGET) query \
		-mode variants \
		-ifn_libraries $(TEST_OUTPUT_DIR)/libraries.tsv \
		-ifn_intervals $(TEST_INTERVALS_SMALL) \
		-ofn_prefix $(TEST_OUTPUT_DIR)/query_variants \
		-min_variants_variant_support 1 \
		-min_variants_library_support 1 \
		-min_variants_coverage_support 1
	@echo "QUERY VARIANTS completed successfully"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

test_coverage: $(TARGET)
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running COVERAGE TEST"
	$(TARGET) coverage \
		-ifn $(TEST_OUTPUT_DIR)/test.aln \
		-ofn_prefix $(TEST_OUTPUT_DIR)/coverage_test
	@echo "COVERAGE TEST completed successfully"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

test_query_all: test_query_full test_query_bin test_query_pileup test_query_consensus test_query_variants test_coverage

########################################################################################
# Test R interface
########################################################################################

# PAF file with dense alignments to test query_full
TEST_DENSE_PAF = examples/align_100_dense.paf
TEST_DENSE_INTERVALS = examples/intervals_dense.txt

# create PAF file with 100 simulated alignments on a 1kb contig
# called once during development
test_create_dense_paf:
	perl pl/simulate_paf.pl $(TEST_DENSE_PAF)
	@echo "PAF file created successfully"

test_R_commands:
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running R SCRIPT TEST"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	Rscript R/example.R \
		$(TEST_DENSE_PAF) \
		$(TEST_DENSE_INTERVALS) \
		$(TEST_BIN_SIZE) \
		$(TEST_OUTPUT_DIR)/R
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

test_R_plot:
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running R PLOT TEST"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	Rscript R/plot.r \
		$(TEST_OUTPUT_DIR)/R_alignments.tsv \
		$(TEST_OUTPUT_DIR)/contig_plot.png
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

test_R: test_R_commands test_R_plot

########################################################################################
# Gene annotation tests  
########################################################################################

# Test gene annotation with dense PAF data
test_genes: $(TARGET)
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running GENE ANNOTATION TEST"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	# First construct ALN file with dense data
	rm -rf $(TEST_OUTPUT_DIR) && mkdir -p $(TEST_OUTPUT_DIR)
	$(TARGET) construct \
		-ifn_paf $(TEST_DENSE_PAF) \
		-ofn $(TEST_OUTPUT_DIR)/test_dense.aln
	@echo "ALN file created, now testing gene annotation..."
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	# Test variants with gene annotation
	@echo "lib_id	aln_fn" > $(TEST_OUTPUT_DIR)/libraries_genes.tsv
	@echo "lib1	$(TEST_OUTPUT_DIR)/test_dense.aln" >> $(TEST_OUTPUT_DIR)/libraries_genes.tsv
	$(TARGET) query \
		-mode variants \
		-ifn_libraries $(TEST_OUTPUT_DIR)/libraries_genes.tsv \
		-ifn_intervals $(TEST_DENSE_INTERVALS) \
		-ofn_prefix $(TEST_OUTPUT_DIR)/variants_genes \
		-ifn_gene_table examples/test_genes.txt \
		-ifn_codon_table table11 \
		-use_genes T \
		-min_variants_variant_support 1 \
		-min_variants_library_support 1 \
		-min_variants_coverage_support 1
	@echo "Variants with gene annotation completed"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	# Test pileup with gene annotation
	$(TARGET) query \
		-ifn_aln $(TEST_OUTPUT_DIR)/test_dense.aln \
		-ifn_intervals $(TEST_DENSE_INTERVALS) \
		-ofn_prefix $(TEST_OUTPUT_DIR)/pileup_genes \
		-mode pileup \
		-ifn_gene_table examples/test_genes.txt \
		-ifn_codon_table table11 \
		-use_genes T \
		-pileup_mode mutated
	@echo "Pileup with gene annotation completed"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	# Test full with gene annotation  
	$(TARGET) query \
		-ifn_aln $(TEST_OUTPUT_DIR)/test_dense.aln \
		-ifn_intervals $(TEST_DENSE_INTERVALS) \
		-ofn_prefix $(TEST_OUTPUT_DIR)/full_genes \
		-mode full \
		-ifn_gene_table examples/test_genes.txt \
		-ifn_codon_table table11 \
		-use_genes T
	@echo "Full query with gene annotation completed"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "Checking output files..."
	@ls -la $(TEST_OUTPUT_DIR)/*genes*
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "GENE ANNOTATION TEST completed successfully"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

########################################################################################
# Test rearrangement detection
########################################################################################

test_rearrange_basic: $(TARGET)
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running REARRANGEMENT BASIC TEST (resolve_seams=no)"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	# First construct ALN file
	rm -rf $(TEST_OUTPUT_DIR) && mkdir -p $(TEST_OUTPUT_DIR)
	$(TARGET) construct \
		-ifn_paf $(TEST_PAF) \
		-ofn $(TEST_OUTPUT_DIR)/test.aln
	@echo "ALN file created, now testing rearrangement detection..."
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	# Test basic rearrangement detection without sequence resolution
	$(TARGET) rearrange \
		-ifn_aln $(TEST_OUTPUT_DIR)/test.aln \
		-ifn_intervals $(TEST_INTERVALS_SMALL) \
		-ofn_prefix $(TEST_OUTPUT_DIR)/rearrange_basic \
		-max_margin 10 \
		-min_element_length 50 \
		-min_anchor_length 200 \
		-max_anchor_mutations_percent 0.01 \
		-max_element_mutation_percent 0.01 \
		-should_verify F \
		-resolve_seams no
	@echo "Checking output files..."
	@ls -la $(TEST_OUTPUT_DIR)/rearrange_basic*
	@echo "Checking CSV structure (first few lines)..."
	@head -5 $(TEST_OUTPUT_DIR)/rearrange_basic.csv || echo "CSV file not found"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "REARRANGEMENT BASIC TEST completed successfully"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

test_rearrange_reference: $(TARGET)
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running REARRANGEMENT REFERENCE TEST (resolve_seams=reference_only)"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	# Construct ALN file with verification enabled
	rm -rf $(TEST_OUTPUT_DIR) && mkdir -p $(TEST_OUTPUT_DIR)
	$(TARGET) construct \
		-ifn_paf $(TEST_PAF) \
		-ifn_reads $(TEST_READS) \
		-ifn_contigs $(TEST_CONTIGS) \
		-verify T \
		-ofn $(TEST_OUTPUT_DIR)/test_verified.aln
	@echo "Verified ALN file created, now testing rearrangement with reference seams..."
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	# Test rearrangement detection with reference sequence resolution
	$(TARGET) rearrange \
		-ifn_aln $(TEST_OUTPUT_DIR)/test_verified.aln \
		-ifn_intervals $(TEST_INTERVALS_SMALL) \
		-ofn_prefix $(TEST_OUTPUT_DIR)/rearrange_reference \
		-max_margin 10 \
		-min_element_length 50 \
		-min_anchor_length 200 \
		-max_anchor_mutations_percent 0.01 \
		-max_element_mutation_percent 0.01 \
		-should_verify T \
		-ifn_contigs $(TEST_CONTIGS) \
		-resolve_seams reference_only
	@echo "Checking output files..."
	@ls -la $(TEST_OUTPUT_DIR)/rearrange_reference*
	@echo "Checking CSV structure (first few lines)..."
	@head -5 $(TEST_OUTPUT_DIR)/rearrange_reference.csv || echo "CSV file not found"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "REARRANGEMENT REFERENCE TEST completed successfully"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

# Combined rearrangement tests
test_rearrange: test_rearrange_basic test_rearrange_reference
	@echo "all rearrangement tests completed successfully"

########################################################################################
# Clean test outputs
########################################################################################

# Clean test outputs
test_clean:
	rm -rf $(TEST_OUTPUT_DIR) 

########################################################################################
# combo rules
########################################################################################

test: test_basic test_full test_query_full test_query_all
	@echo "all tests completed successfully"
