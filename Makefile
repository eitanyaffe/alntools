CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra
LDFLAGS = -lz

# Detect OS first
UNAME_S := $(shell uname -s)

# Metal support for Apple Silicon GPUs (enabled when Xcode tools + Metal GPU available)
METAL_SUPPORT := false
ifeq ($(UNAME_S),Darwin)
    # Check both Xcode tools AND actual Metal GPU support
    XCODE_AVAILABLE := $(shell xcode-select -p >/dev/null 2>&1 && echo true || echo false)
    METAL_GPU_AVAILABLE := $(shell system_profiler SPDisplaysDataType | grep -q "Metal" && echo true || echo false)
    METAL_AVAILABLE := $(shell [ "$(XCODE_AVAILABLE)" = "true" ] && [ "$(METAL_GPU_AVAILABLE)" = "true" ] && echo true || echo false)
    ifeq ($(METAL_AVAILABLE),true)
        METAL_SUPPORT := true
        CXXFLAGS += -DMETAL_SUPPORT
        LDFLAGS += -framework Metal -framework MetalPerformanceShaders -framework Foundation
        OBJCXXFLAGS = -std=c++17 -Wall -Wextra -DMETAL_SUPPORT -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
    endif
endif

# OpenMP support - different for different systems
ifeq ($(UNAME_S),Darwin)
    # macOS - check if OpenMP is available
    OPENMP_TEST := $(shell echo | $(CXX) -fopenmp -x c++ -c - -o /dev/null 2>/dev/null && echo true || echo false)
    ifeq ($(OPENMP_TEST),true)
        CXXFLAGS += -fopenmp
        LDFLAGS += -fopenmp
    else
        # Try with libomp if available (Homebrew installation)
        HOMEBREW_LIBOMP := $(shell test -d /opt/homebrew/opt/libomp && echo true || echo false)
        ifeq ($(HOMEBREW_LIBOMP),true)
            CXXFLAGS += -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
            LDFLAGS += -L/opt/homebrew/opt/libomp/lib -lomp
        else
            # Try with system libomp
            LIBOMP_TEST := $(shell echo | $(CXX) -Xpreprocessor -fopenmp -lomp -x c++ -c - -o /dev/null 2>/dev/null && echo true || echo false)
            ifeq ($(LIBOMP_TEST),true)
                CXXFLAGS += -Xpreprocessor -fopenmp
                LDFLAGS += -lomp
            else
                $(warning OpenMP not available, compiling without threading support)
            endif
        endif
    endif
else
    # Linux - assume OpenMP is available
    CXXFLAGS += -fopenmp
    LDFLAGS += -fopenmp
endif

SRC_DIR = cpp

ifeq ($(UNAME_S),Linux)
    OBJ_DIR = obj/linux
    BIN_DIR = bin/linux
else ifeq ($(UNAME_S),Darwin)
    OBJ_DIR = obj/macos
    BIN_DIR = bin/macos
endif

# Header files
HDRS = $(wildcard $(SRC_DIR)/*.h)

# Source files
SRC_FILES_TO_EXCLUDE := $(SRC_DIR)/aln_R.cpp
SRCS = $(filter-out $(SRC_FILES_TO_EXCLUDE), $(wildcard $(SRC_DIR)/*.cpp))

# Metal Objective-C++ files
ifeq ($(METAL_SUPPORT),true)
    METAL_SRCS = $(SRC_DIR)/QueryBin_metal.mm
    METAL_OBJS = $(METAL_SRCS:$(SRC_DIR)/%.mm=$(OBJ_DIR)/%.o)
else
    METAL_SRCS =
    METAL_OBJS =
endif
# Define OBJS based on the filtered SRCS list
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o) $(METAL_OBJS)

# Binary name
TARGET = $(BIN_DIR)/alntools

# Installation directory - local user directory
INSTALL_DIR = $(HOME)/.local/bin

# Default target
all: $(TARGET)

# Create necessary directories
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Compile source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HDRS) | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile Metal Objective-C++ files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.mm $(HDRS) | $(OBJ_DIR)
	$(CXX) $(OBJCXXFLAGS) -c $< -o $@

# Link object files
$(TARGET): $(OBJS) | $(BIN_DIR)
	$(CXX) $(OBJS) -o $(TARGET) $(LDFLAGS)

# Install target
install: $(TARGET)
	mkdir -p $(INSTALL_DIR)
	cp $(TARGET) $(INSTALL_DIR)/

# Include test rules
include test.mk

# Clean build files
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# Phony targets
.PHONY: all clean install 
