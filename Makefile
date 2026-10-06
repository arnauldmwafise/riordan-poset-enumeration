# Makefile -- builds everything without CMake. `make` builds, `make test` runs all tests,
# `make check-rcount` checks the large-scale counter against verified values.
CXX      ?= g++
CC       ?= gcc
CXXFLAGS ?= -std=c++20 -O2 -Wall -Wextra -Wpedantic
INC       = -Iinclude
HDRS      = $(wildcard include/riordan/*.hpp)

# nauty (third_party/nauty, Apache-2.0): built for graphs with at most 64 vertices,
# one 64-bit word per row, with thread-local storage so rcount can use threads.
NAUTY_DIR   = third_party/nauty
NAUTY_FLAGS = -DMAXN=WORDSIZE -DWORDSIZE=64
NAUTY_SRCS  = nauty nautil naugraph schreier naurng nautinv
NAUTY_OBJS  = $(addprefix bin/nauty/,$(addsuffix .o,$(NAUTY_SRCS)))

TEST_SRCS = $(wildcard tests/test_*.cpp)
TEST_BINS = $(patsubst tests/%.cpp,bin/tests/%,$(TEST_SRCS))

all: bin/riordan bin/rcount bin/independent_rcount bin/quickstart bin/sample_listing $(TEST_BINS)

bin/riordan: tools/riordan_cli.cpp $(HDRS) | bin
	$(CXX) $(CXXFLAGS) $(INC) $< -o $@
bin/quickstart: examples/quickstart.cpp $(HDRS) | bin
	$(CXX) $(CXXFLAGS) $(INC) $< -o $@
bin/sample_listing: examples/sample_listing.cpp $(HDRS) | bin
	$(CXX) $(CXXFLAGS) $(INC) $< -o $@
bin/tests/%: tests/%.cpp tests/test_common.hpp $(HDRS) | bin
	mkdir -p bin/tests
	$(CXX) $(CXXFLAGS) $(INC) $< -o $@

bin/nauty/%.o: $(NAUTY_DIR)/%.c | bin
	mkdir -p bin/nauty
	$(CC) -O3 $(NAUTY_FLAGS) -c $< -o $@
bin/rcount: tools/rcount.cpp $(HDRS) $(NAUTY_OBJS) | bin
	$(CXX) -std=c++20 -O3 -Wall $(NAUTY_FLAGS) $(INC) -I$(NAUTY_DIR) $< $(NAUTY_OBJS) -o $@ -pthread

# Independent verification counter: compiled ONLY against verify/independent (not the library).
bin/independent_rcount: verify/independent_rcount.cpp $(wildcard verify/independent/include/riordan/*.hpp) $(NAUTY_OBJS) | bin
	$(CXX) -std=c++20 -O3 -Wall -Wextra $(NAUTY_FLAGS) -Iverify/independent/include -I$(NAUTY_DIR) $< $(NAUTY_OBJS) -o $@ -pthread

bin:
	mkdir -p bin

test: $(TEST_BINS)
	@fail=0; for t in $(TEST_BINS); do echo "== $$t"; $$t | tail -1 || fail=1; $$t > /dev/null || fail=1; done; \
	if [ $$fail = 0 ]; then echo "all test suites passed"; else echo "SOME TEST SUITES FAILED"; exit 1; fi

# Compare rcount (both canonical forms, and --passes mode) with the verified values r(1..16).
check-rcount: bin/rcount bin/independent_rcount
	@sh scripts/check_rcount.sh ./bin/rcount
	@sh scripts/check_independent.sh ./bin/independent_rcount

clean:
	rm -rf bin

.PHONY: all test clean check-rcount
