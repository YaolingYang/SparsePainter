# Makefile for SparsePainter

# Compiler and compiler flags
CXX = g++
CXXFLAGS = -std=c++0x -g -O3

# Libraries
LIBS = -lz -fopenmp -larmadillo

# Source file and target executable
SRC = SparsePainter.cpp
TARGET = SparsePainter

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LIBS)

clean:
	rm -f $(TARGET)

.PHONY: all clean
