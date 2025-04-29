# Test rules for alntools
# This file contains test rules that are included in the main Makefile

# input
TEST_PAF = examples/align_100.paf
TEST_READS = examples/reads_100.fq
TEST_CONTIGS = examples/contigs_100.fa

# output
TEST_BASIC_OUTPUT_DIR = output/basic
TEST_FULL_OUTPUT_DIR = output/full

TEST_EXAMPLE_OUTPUT_DIR = output/example

TEST_UNCOMPRESSED = examples/.uncompressed

TEST_INTERVALS_SMALL = examples/intervals_small.txt
TEST_INTERVALS_LARGE = examples/intervals_large.txt

# Test targets
.PHONY: test test_basic test_full test_query clean-test

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
	rm -rf $(TEST_FULL_OUTPUT_DIR) && mkdir -p $(TEST_FULL_OUTPUT_DIR)
	@echo "running FULL TEST, with paf validation"
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

# construct ALN with validation
create_example: $(TARGET)
	rm -rf $(TEST_EXAMPLE_OUTPUT_DIR) && mkdir -p $(TEST_EXAMPLE_OUTPUT_DIR)
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running EXAMPLE TEST"
	$(TARGET) verify \
		-ifn_aln $(TEST_BASIC_OUTPUT_DIR)/test.aln \
		-ifn_reads $(TEST_READS) \
		-ifn_contigs $(TEST_CONTIGS) \
		-max_reads 100 \
		-ofn_reads examples/short_reads.fq \
		-ofn_contigs examples/short_contigs.fa

test_query_full: $(TARGET)
	$(TARGET) query \
		-ifn_aln $(TEST_BASIC_OUTPUT_DIR)/test.aln \
		-ifn_intervals $(TEST_INTERVALS_LARGE) \
		-ofn_prefix $(TEST_BASIC_OUTPUT_DIR)/query \
		-mode full

test_query_bin: $(TARGET)
	$(TARGET) query \
		-ifn_aln $(TEST_BASIC_OUTPUT_DIR)/test.aln \
		-ifn_intervals $(TEST_INTERVALS_LARGE) \
		-ofn_prefix $(TEST_BASIC_OUTPUT_DIR)/query \
		-mode bin \
		-binsize 10000 \
		-skip_empty_bins F

test_query_pileup: $(TARGET)
	$(TARGET) query \
		-ifn_aln $(TEST_BASIC_OUTPUT_DIR)/test.aln \
		-ifn_intervals $(TEST_INTERVALS_SMALL) \
		-ofn_prefix $(TEST_BASIC_OUTPUT_DIR)/query \
		-mode pileup \
		-pileup_mode mutated

# Run all tests
test: test_basic test_full test_query
	@echo "all tests completed successfully"

# Clean test outputs
test_clean:
	rm -rf $(TEST_BASIC_OUTPUT_DIR) $(TEST_FULL_OUTPUT_DIR) 
