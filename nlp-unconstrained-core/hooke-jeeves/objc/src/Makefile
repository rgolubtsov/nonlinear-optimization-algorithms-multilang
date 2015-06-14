# =============================================================================
# nlp-unconstrained-core/hooke-jeeves/objc/src/Makefile
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
TMPS    = *.d
BIN_DIR = ../bin
OBJ_DIR = ../obj

# Specify flags and other vars here.
CFLAGS = -Wall `gnustep-config --objc-flags`
LDLIBS = -lobjc -lgnustep-base

MV      = mv
MVFLAGS = -v

# Making the first target (rosenbrock).
$(DEPS1) $(DEPSN): CFLAGS += -c -o
$(DEPS1) $(DEPSN): %.o: %.m
	$(CC) $(CFLAGS) $@ $<

$(EXEC1): $(DEPS1) $(DEPSN)

# Making the second target (woods).
$(DEPS2): CFLAGS += -DWOODS -c -o
$(DEPS2): %$(EXEC2S).o: %.m
	$(CC) $(CFLAGS) $@ $<

$(EXEC2): $(DEPS2) $(DEPSN)

.PHONY: rosenbrock woods all clean

rosenbrock: $(EXEC1)

woods: $(EXEC2)

all: rosenbrock woods
	$(MV) $(MVFLAGS) $(EXEC1) $(EXEC2) $(BIN_DIR)
	$(MV) $(MVFLAGS) $(DEPS1) $(DEPS2) $(DEPSN) $(TMPS) $(OBJ_DIR)

clean:
	$(RM) $(EXEC1) $(DEPS1) $(EXEC2) $(DEPS2) $(DEPSN) $(TMPS)

# =============================================================================
# vim:set nu:et:ts=4:sw=4:
# =============================================================================
