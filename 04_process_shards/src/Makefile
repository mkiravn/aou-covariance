# Makefile for grm_class_tool
#
#   make        -- build grm_class_tool
#   make test   -- build, then run the correctness test against a Python reference
#   make clean  -- remove build + test artifacts

CXX      = g++
CXXFLAGS = -O2 -Wall -Wextra -std=c++17
BIN      = grm_class_tool

.PHONY: all test clean

all: $(BIN)

$(BIN): grm_class_tool.cpp
	$(CXX) $(CXXFLAGS) grm_class_tool.cpp -o $(BIN)

test: $(BIN)
	./test/run_test.sh

clean:
	rm -f $(BIN)
	rm -rf test/tmp
