# Compiler
CXX = g++

# Directories and Libraries
ARMA_DIR = ./armadillo
INCLUDE_DIR = $(ARMA_DIR)/include
LIB_DIR = $(ARMA_DIR)

# Compiler Flags
CXXFLAGS = -I$(INCLUDE_DIR) -std=c++0x -g -O3 -fopenmp
LDFLAGS = -L$(LIB_DIR)

# Determine OS type and set linker flags accordingly
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    LDLIBS = -larmadillo -lblas -llapack -lz -lpthread -Wl,-rpath=$(LIB_DIR)
endif
ifeq ($(UNAME_S),Darwin)  # Darwin is the result for macOS
    ## We will set CXX and LDLIBS specially
    $(info CXX = $(CXX))
    $(info On Mac, we cannot use the Xcode compiler based on clang. Checking for the correct version...)
    CXXVERSION=$($(CXX) -v 2>&1 | grep gcc) # check for gcc c++
    ifeq (,$(findstring gcc,$(CXXVERSION))) # Not found gcc
        $(warning The specified CXX is not supported. Searching for a gcc version in your PATH...) 
        CXX=$(shell ./whereisglob g++*) # search for g++
        $(info Found CXX = $(CXX))
    endif
    CXXVERSION=$($(CXX) -v 2>&1 | grep gcc)
    ifneq (,$(findstring gcc,$(CXXVERSION))) # Not found gcc
        $(error You don't have a gcc supported compiler as the first g++* in your path. You must a. have installed g++, and b. specify the version my running make CXX=g++-<version>, for example, g++-13)
    endif
    $(info Proceeding with gcc g++ executable $(CXX))
    LDLIBS = -larmadillo -framework Accelerate -lz -lpthread -Wl,-ld_classic -Wl,-rpath $(LIB_DIR)
endif

# Target
TARGET = SparsePainter

all: $(TARGET)

$(TARGET): SparsePainter.cpp
	$(CXX) $< -o $@ $(CXXFLAGS) $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(TARGET)
