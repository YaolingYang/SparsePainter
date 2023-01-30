## Compile all pbwt and dpbwt files

CC = g++
CPPFLAGS = -O3 -std=c++17


## Detect if we're on windows
ifeq ($(OS),Windows_NT) 
    detected_OS := Windows
else
    detected_OS := $(shell sh -c 'uname 2>/dev/null || echo Unknown')
endif

objects : lmq_3swp_pbwt lmq_3swp_dpbwt lmq_1swp_dpbwt lmq_1swp_pbwt

all: $(objects)

$(objects): %: %.cpp
	$(CC) $(CFLAGS) -o $@ $<
