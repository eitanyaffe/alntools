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

.PHONY: test test_basic test_full test_query_full test_query_bin \
test_query_pileup test_query_all test_R_all test_R_commands test_R_plot \
test_create_dense_paf clean-test test-r-load

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

test_query_all: test_query_full test_query_bin test_query_pileup

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

test_R_all: test_R_commands test_R_plot
########################################################################################
# combo rules
########################################################################################

test: test_basic test_full test_query_full test_query_all
	@echo "all tests completed successfully"

# Clean test outputs
test_clean:
	rm -rf $(TEST_OUTPUT_DIR) 
