# Test rules for alntools
# This file contains test rules that are included in the main Makefile

# input
TEST_PAF = examples/short_align.paf
TEST_READS = examples/short_reads.fq
TEST_CONTIGS = examples/short_contigs.fa

# output
TEST_BASIC_OUTPUT_DIR = output/basic
TEST_FULL_OUTPUT_DIR = output/full

TEST_EXAMPLE_OUTPUT_DIR = output/example

TEST_UNCOMPRESSED = examples/.uncompressed

# Test targets
.PHONY: test test_basic test_full clean-test

# construct ALN w/o validation
test_basic: $(TARGET)
	rm -rf $(TEST_BASIC_OUTPUT_DIR) && mkdir -p $(TEST_BASIC_OUTPUT_DIR)
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	@echo "running BASIC TEST, without paf validation"
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
	rm -rf $(TEST_FULL_OUTPUT_DIR) && mkdir -p $(TEST_FULL_OUTPUT_DIR)
	@echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
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

# Run all tests
test: test_basic test_full
	@echo "all tests completed successfully"

# Clean test outputs
test_clean:
	rm -rf $(TEST_BASIC_OUTPUT_DIR) $(TEST_FULL_OUTPUT_DIR) 
