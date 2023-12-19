# Compiler settings - Can be customized.
CXX = g++
CXXFLAGS = -Wall -std=c++11 -pg -Iinclude -Ilib -O3 -fopenmp 

# Build directories
SRCDIR = src
BUILDDIR = build
TESTDIR = tests

# Target executable for tests
TESTTARGET = $(BUILDDIR)/testInputData 

# Source and Object files
SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(patsubst $(SRCDIR)/%.cpp, $(BUILDDIR)/%.o, $(SRCS))
TESTSRCS = $(wildcard $(TESTDIR)/*.cpp)
TESTOBJS = $(filter-out $(BUILDDIR)/main.o, $(OBJS)) $(patsubst $(TESTDIR)/%.cpp, $(BUILDDIR)/%.o, $(TESTSRCS))

# Default rule
all: test

# Test build rule
test: $(TESTTARGET)
	@echo "Running tests..."
	@./$(TESTTARGET) 

$(TESTTARGET): $(TESTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# General rule for building objects
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILDDIR)/%.o: $(TESTDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILDDIR):
	@mkdir -p $(BUILDDIR)

# Clean up
clean:
	rm -rf $(BUILDDIR)/* $(TESTTARGET)

.PHONY: all test clean

