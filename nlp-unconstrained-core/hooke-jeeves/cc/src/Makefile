# =============================================================================
# nlp-unconstrained-core/hooke-jeeves/cc/src/Makefile
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
# =============================================================================

EXEC1   = hooke
EXEC2S  = -woods
EXEC2   = $(EXEC1)$(EXEC2S)
DEPS1S  = rosenbrock
DEPS1   = $(EXEC1).o $(DEPS1S).o
DEPS2S  = woods
DEPS2   = $(EXEC2).o $(DEPS2S)$(EXEC2S).o
DEPSN1  = funevals
DEPSN   = $(DEPSN1).o
BIN_DIR = ../bin
OBJ_DIR = ../obj

# Specify flags and other vars here.
# The following standards-compliance options work well with g++ 5.1.0.
#CXX_STD = c++03 # Same as c++98, same as -ansi .
#CXX_STD = c++11
#CXX_STD = c++14

# Travis CI g++ (4.6.3) doesn't know about the -std=c++14 (and c++11) option.
# Using the following option instead.
CXX_STD = c++0x

CXXFLAGS = -Wall -pedantic -std=$(CXX_STD) -O2
LDLIBS   = -lstdc++

MV      = mv
MVFLAGS = -v

# Making the first target (rosenbrock).
$(DEPS1) $(DEPSN): %.o: %.cc

$(EXEC1): $(DEPS1) $(DEPSN)

# Making the second target (woods).
$(DEPS2): CXXFLAGS += -DWOODS -c -o
$(DEPS2): %$(EXEC2S).o: %.cc
	$(CXX) $(CXXFLAGS) $@ $<

$(EXEC2): $(DEPS2) $(DEPSN)

.PHONY: rosenbrock woods all clean

rosenbrock: $(EXEC1)

woods: $(EXEC2)

all: rosenbrock woods
	$(MV) $(MVFLAGS) $(EXEC1) $(EXEC2) $(BIN_DIR)
	$(MV) $(MVFLAGS) $(DEPS1) $(DEPS2) $(DEPSN) $(OBJ_DIR)

clean:
	$(RM) $(EXEC1) $(DEPS1) $(EXEC2) $(DEPS2) $(DEPSN)

# =============================================================================
# vim:set nu:et:ts=4:sw=4:
# =============================================================================
