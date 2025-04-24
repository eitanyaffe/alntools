CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra
LDFLAGS = -lz

SRC_DIR = cpp

# Detect the operating system
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Linux)
    OBJ_DIR = obj/linux
    BIN_DIR = bin/linux
else ifeq ($(UNAME_S),Darwin)
    OBJ_DIR = obj/macos
    BIN_DIR = bin/macos
endif

# Source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# Binary name
TARGET = $(BIN_DIR)/alntools

# Default target
all: $(TARGET)

# Create necessary directories
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Compile source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link object files
$(TARGET): $(OBJS) | $(BIN_DIR)
	$(CXX) $(OBJS) -o $(TARGET) $(LDFLAGS)

# Include test rules
include mk/test.mk

# Clean build files
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# Phony targets
.PHONY: all clean 
