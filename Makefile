# Compiler
CXX = g++

# Directories and Libraries
ARMA_DIR = ./armadillo
INCLUDE_DIR = $(ARMA_DIR)/include
LIB_DIR = $(ARMA_DIR)

# Compiler Flags
CXXFLAGS = -I$(INCLUDE_DIR) -std=c++0x -g -O3 -fopenmp
LDFLAGS = -L$(LIB_DIR) -Wl,-rpath=$(LIB_DIR)

# Determine OS type and set linker flags accordingly
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    LDLIBS = -larmadillo -lblas -llapack -lz -lpthread
endif
ifeq ($(UNAME_S),Darwin)  # Darwin is the result for macOS
    LDLIBS = -larmadillo -framework Accelerate -lz -lpthread
endif

# Target
TARGET = SparsePainter

all: $(TARGET)

$(TARGET): SparsePainter.cpp
	$(CXX) $< -o $@ $(CXXFLAGS) $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(TARGET)
