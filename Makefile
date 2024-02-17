# Target
TARGET = SparsePainter

all: armadillo $(TARGET)

# Compile armadillo 12.6.5

#.PHONY: setup_armadillo

armadillo:
	tar -xf armadillo-12.6.5.tar.xz && rm armadillo-12.6.5.tar.xz

# Compiler
CXX = g++

# Directories and Libraries
ARMA_DIR = ./armadillo-12.6.5
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
    $(shell chmod +x ./whereisglob)
    CXXVERSION=$(shell $(CXX) -v 2>&1 | grep Apple) # check for Apple c++
    ifeq (Apple,$(findstring Apple,$(CXXVERSION))) # Found apple c++
        $(info On Mac, we cannot use the Xcode compiler based on clang (called Darwin).)
        $(warning The specified CXX is not supported. Searching for another version in your PATH...) 
        CXX=$(shell ./whereisglob g++*) # search for g++
        $(info Found CXX = $(CXX))
    endif
    CXXVERSION=$(shell $(CXX) -v 2>&1 | grep Apple) # check for Apple c++
    ifeq (Apple,$(findstring Apple,$(CXXVERSION))) # Found Apple c++
        $(error You have a native Apple Xcode compiler as the first g++* in your path. A gcc compiler is strongly recommended. a. have installed g++, and b. specify the version my running make CXX=g++-<version>, for example, CXX=g++-13)
    endif
    $(info Proceeding with gcc g++ executable $(CXX))
    LDLIBS = -larmadillo -framework Accelerate -lz -lpthread -Wl,-ld_classic -Wl,-rpath $(LIB_DIR)
endif

$(TARGET): SparsePainter.cpp
	$(CXX) $< -o $@ $(CXXFLAGS) $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(TARGET)
