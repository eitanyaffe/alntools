# Test rules for alntools
# This file contains test rules that can be included in the main Makefile

# Directories
EXAMPLES_DIR = examples
OUTPUT_DIR = output

# Test targets
.PHONY: test test-construct test-info test-save clean-test

# Create output directory
$(OUTPUT_DIR):
	mkdir -p $(OUTPUT_DIR)

# Test construct command
test-construct: $(TARGET) $(OUTPUT_DIR)
	@echo "Testing PAF file construction..."
	@$(TARGET) construct $(EXAMPLES_DIR)/example.paf $(OUTPUT_DIR)/example.aln
	@echo "Test completed successfully!"

# Test info command
test-info: $(TARGET) $(OUTPUT_DIR)
	@echo "Testing info command..."
	@$(TARGET) info $(OUTPUT_DIR)/example.aln
	@echo "Test completed successfully!"

# Test save command
test-save: $(TARGET) $(OUTPUT_DIR)
	@echo "Testing save command..."
	@$(TARGET) save $(OUTPUT_DIR)/example.aln $(OUTPUT_DIR)/example
	@echo "Test completed successfully!"

# Run all tests
test: test-construct test-info test-save
	@echo "All tests completed successfully!"

# Clean test outputs
clean-test:
	rm -rf $(OUTPUT_DIR) 