# =============================================================================
# nlp-unconstrained-core/hooke-jeeves/objc/Makefile
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
# =============================================================================
# Copyright (C) 2015 Radislav (Radic) Golubtsov

EXEC1  = hooke
EXEC2S = -woods
EXEC2  = $(EXEC1)$(EXEC2S)
DEPS1S = rosenbrock
DEPS1  = $(EXEC1).o $(DEPS1S).o
DEPS2S = woods
DEPS2  = $(EXEC2).o $(DEPS2S)$(EXEC2S).o
TMPS   = *.d

# Specify flags and other vars here.
CFLAGS = -Wall `gnustep-config --objc-flags`
LDLIBS = -lobjc -lgnustep-base

# Making the first target (rosenbrock).
$(DEPS1): CFLAGS += -c -o
$(DEPS1): %.o: %.m
	$(CC) $(CFLAGS) $@ $<

$(EXEC1): $(DEPS1)

# Making the second target (woods).
$(DEPS2): CFLAGS += -DWoods -c -o
$(DEPS2): %$(EXEC2S).o: %.m
	$(CC) $(CFLAGS) $@ $<

$(EXEC2): $(DEPS2)

.PHONY: rosenbrock woods all clean

rosenbrock: $(EXEC1)

woods: $(EXEC2)

all: rosenbrock woods

clean:
	$(RM) $(EXEC1) $(DEPS1) $(EXEC2) $(DEPS2) $(TMPS)

# =============================================================================
# vim:set nu:et:ts=4:sw=4:
# =============================================================================
