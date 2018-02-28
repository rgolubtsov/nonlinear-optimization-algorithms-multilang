#
# nlp-unconstrained-cli/nelder-mead/__orig/src/Makefile
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Nelder-Mead nonlinear unconstrained minimization algorithm.
# =============================================================================
# Written by Radislav (Radicchio) Golubtsov, 2015-2017
#
# This is free and unencumbered software released into the public domain.
#
# Anyone is free to copy, modify, publish, use, compile, sell, or
# distribute this software, either in source code form or as a compiled
# binary, for any purpose, commercial or non-commercial, and by any
# means.
#
# (See the LICENSE file at the top of the source tree.)
#

EXEC1   = nelmin
EXEC2S  = -woods
EXEC2   = $(EXEC1)$(EXEC2S)
DEPS1   = $(EXEC1).o
DEPS2   = $(EXEC2).o
BIN_DIR = ../bin
OBJ_DIR = ../obj

# Specify flags and other vars here.
CC      = gfortran
CFLAGS  = -Wall -cpp -ffixed-form -pedantic -std=legacy -O3

MV      = mv
MVFLAGS = -v

# Making the first target (rosenbrock).
$(DEPS1): CFLAGS += -c -o
$(DEPS1): %.o: %.f
	$(CC) $(CFLAGS) $@ $<

$(EXEC1): $(DEPS1)

# Making the second target (woods).
$(DEPS2): CFLAGS += -DWOODS -c -o
$(DEPS2): %$(EXEC2S).o: %.f
	$(CC) $(CFLAGS) $@ $<

$(EXEC2): $(DEPS2)

.PHONY: rosenbrock woods all clean

rosenbrock: $(EXEC1)

woods: $(EXEC2)

all: rosenbrock woods
	$(MV) $(MVFLAGS) $(EXEC1) $(EXEC2) $(BIN_DIR)
	$(MV) $(MVFLAGS) $(DEPS1) $(DEPS2) $(OBJ_DIR)

clean:
	$(RM) $(EXEC1) $(DEPS1) $(EXEC2) $(DEPS2)

# vim:set nu ts=4 sw=4:
