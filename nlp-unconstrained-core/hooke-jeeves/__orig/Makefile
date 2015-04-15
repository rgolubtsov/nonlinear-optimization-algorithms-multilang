# =============================================================================
# nlp-unconstrained-core/hooke-jeeves/__orig/Makefile
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 1.0
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
DEPS1  = $(EXEC1).o
DEPS2  = $(EXEC2).o

# == Specify flags and other vars here.
# To use gcc instead of tcc, simply comment out next line, just like this:
# "#CC = tcc" - and GNU C Compiler will be used automatically.
CC     = tcc
CFLAGS = -Wall
LDLIBS = -lm

# == Making the first target (rosenbrock).
$(DEPS1): %.o: %.c

$(EXEC1): $(DEPS1)

# == Making the second target (woods).
$(DEPS2): CFLAGS += -DWoods -c -o
$(DEPS2): %$(EXEC2S).o: %.c
	$(CC) $(CFLAGS) $@ $<

$(EXEC2): $(DEPS2)

.PHONY: rosenbrock woods all clean

rosenbrock: $(EXEC1)

woods: $(EXEC2)

all: rosenbrock woods

clean:
	$(RM) $(EXEC1) $(DEPS1) $(EXEC2) $(DEPS2)

# =============================================================================
# vim:set nu:et:ts=4:sw=4:
# =============================================================================
