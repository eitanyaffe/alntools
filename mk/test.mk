# Test rules for alntools
# This file contains test rules that are included in the main Makefile

# input
TEST_PAF = examples/align_100.paf
TEST_READS = examples/reads_100.fq
TEST_CONTIGS = examples/contigs_100.fa

# output
TEST_BASIC_OUTPUT_DIR = output/basic
TEST_FULL_OUTPUT_DIR = output/full

TEST_INTERVALS_SMALL = examples/intervals_small.txt
TEST_INTERVALS_LARGE = examples/intervals_large.txt

TEST_BIN_SIZE = 1000
# Test targets
.PHONY: test test_basic test_full test_query clean-test test-r-load

########################################################################################
# creating aln files
########################################################################################

# construct ALN w/o validation
test_basic: $(TARGET)
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running BASIC TEST, without paf validation"
	rm -rf $(TEST_BASIC_OUTPUT_DIR) && mkdir -p $(TEST_BASIC_OUTPUT_DIR)
	$(TARGET) construct \
		-ifn_paf $(TEST_PAF) \
		-ofn $(TEST_BASIC_OUTPUT_DIR)/test.aln
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	$(TARGET) info \
		-ifn $(TEST_BASIC_OUTPUT_DIR)/test.aln
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	$(TARGET) extract \
		-ifn $(TEST_BASIC_OUTPUT_DIR)/test.aln \
		-ofn_prefix $(TEST_BASIC_OUTPUT_DIR)/test
	@echo "BASIC TEST completed successfully"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

# construct ALN with validation
test_full: $(TARGET)
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running FULL TEST, with paf validation"
	rm -rf $(TEST_FULL_OUTPUT_DIR) && mkdir -p $(TEST_FULL_OUTPUT_DIR)
	$(TARGET) construct \
		-ifn_paf $(TEST_PAF) \
		-ifn_reads $(TEST_READS) \
		-ifn_contigs $(TEST_CONTIGS) \
		-verify T \
		-ofn $(TEST_FULL_OUTPUT_DIR)/test.aln
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	$(TARGET) info \
		-ifn $(TEST_FULL_OUTPUT_DIR)/test.aln
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	$(TARGET) extract \
		-ifn $(TEST_FULL_OUTPUT_DIR)/test.aln \
		-ofn_prefix $(TEST_FULL_OUTPUT_DIR)/test
	@echo "FULL TEST completed successfully"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

########################################################################################
# query tests
########################################################################################

test_query_full: $(TARGET)
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running QUERY FULL"
	$(TARGET) query \
		-ifn_aln $(TEST_BASIC_OUTPUT_DIR)/test.aln \
		-ifn_intervals $(TEST_INTERVALS_LARGE) \
		-ofn_prefix $(TEST_BASIC_OUTPUT_DIR)/query \
		-mode full
	@echo "QUERY FULL completed successfully"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

test_query_bin: $(TARGET)
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running QUERY BIN"
	$(TARGET) query \
		-ifn_aln $(TEST_BASIC_OUTPUT_DIR)/test.aln \
		-ifn_intervals $(TEST_INTERVALS_SMALL) \
		-ofn_prefix $(TEST_BASIC_OUTPUT_DIR)/query \
		-mode bin \
		-binsize $(TEST_BIN_SIZE)
	@echo "QUERY BIN completed successfully"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

test_query_pileup: $(TARGET)
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running QUERY PILEUP"
	$(TARGET) query \
		-ifn_aln $(TEST_BASIC_OUTPUT_DIR)/test.aln \
		-ifn_intervals $(TEST_INTERVALS_SMALL) \
		-ofn_prefix $(TEST_BASIC_OUTPUT_DIR)/query \
		-mode pileup \
		-pileup_mode mutated
	@echo "QUERY PILEUP completed successfully"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

test_query_all: test_query_full test_query_bin test_query_pileup

########################################################################################
# Demonstrate R interaface
########################################################################################

test_R:
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running R SCRIPT TEST"
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	Rscript R/example.R \
		$(TEST_BASIC_OUTPUT_DIR)/test.aln \
		$(TEST_INTERVALS_SMALL) \
		$(TEST_BIN_SIZE) \
		$(TEST_BASIC_OUTPUT_DIR)/R
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

########################################################################################
# combo rules
########################################################################################

test: test_basic test_full test_query_full test_query_all test_R
	@echo "all tests completed successfully"

# Clean test outputs
test_clean:
	rm -rf $(TEST_BASIC_OUTPUT_DIR) $(TEST_FULL_OUTPUT_DIR) 
