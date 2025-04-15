# Test rules for alntools
# This file contains test rules that can be included in the main Makefile

# Directories
TEST_INPUT_FILE = examples/input.paf
TEST_OUTPUT_DIR = output

# Test targets
.PHONY: test test-construct test-info test-save clean-test

# Create output directory
$(TEST_OUTPUT_DIR):

# Test construct command
test-construct: $(TARGET) $(TEST_OUTPUT_DIR)
	@echo "Testing PAF file construction..."
	$(TARGET) construct \
		-ifn $(TEST_INPUT_FILE) \
		-ofn $(TEST_OUTPUT_DIR)/example.aln
	@echo "Test completed successfully!"

# Test info command
test-info: $(TARGET) $(TEST_OUTPUT_DIR)
	@echo "Testing info command..."
	$(TARGET) info \
		-ifn $(TEST_OUTPUT_DIR)/example.aln
	@echo "Test completed successfully!"

# Test save command
test-save: $(TARGET) $(TEST_OUTPUT_DIR)
	@echo "Testing save command..."
	$(TARGET) save \
		-ifn $(TEST_OUTPUT_DIR)/example.aln \
		-ofn_prefix $(TEST_OUTPUT_DIR)/example
	@echo "Test completed successfully!"

# Run all tests
test: test-construct test-info test-save
	@echo "All tests completed successfully!"

# Clean test outputs
clean-test:
	rm -rf $(TEST_OUTPUT_DIR) 
