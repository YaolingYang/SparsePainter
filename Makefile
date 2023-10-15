# Compiler
CXX = g++

# Directories and Libraries
ARMA_DIR = ./armadillo
INCLUDE_DIR = $(ARMA_DIR)/include
LIB_DIR = $(ARMA_DIR)

# Compiler and Linker Flags
CXXFLAGS = -I$(INCLUDE_DIR) -std=c++0x -g -O3 -fopenmp
LDFLAGS = -L$(LIB_DIR) -Wl,-rpath=$(LIB_DIR)
LDLIBS = -larmadillo -llapack -lblas -lz -lpthread

# Target
TARGET = SparsePainter

all: $(TARGET)

$(TARGET): SparsePainter.cpp
	$(CXX) $< -o $@ $(CXXFLAGS) $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(TARGET)
