#
# nlp-unconstrained-cli/hooke-jeeves/vala/src/Makefile
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
#
# (See outer Makefile to understand how this one is processed.)
# =============================================================================
# Written by Radislav (Radicchio) Golubtsov, 2016
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

EXEC1   = hooke
EXEC2   = $(EXEC1)-woods
SRCS1   = $(EXEC1).vala   \
          rosenbrock.vala \
          $(SRCSN)
SRCS2   = $(EXEC1).vala   \
          woods.vala      \
          $(SRCSN)
SRCSN   = funevals.vala
BIN_DIR = ../bin

# Specify flags and other vars here.
VC      = valac
VCFLAGS = --target-glib=2.40 --pkg=posix

MV      = mv
MVFLAGS = -v

# Making the first target (rosenbrock).
$(EXEC1): VCFLAGS += -o
$(EXEC1): $(SRCS1)
	$(VC) $(VCFLAGS) $(EXEC1) $(SRCS1)

# Making the second target (woods).
$(EXEC2): VCFLAGS += -D WOODS -o
$(EXEC2): $(SRCS2)
	$(VC) $(VCFLAGS) $(EXEC2) $(SRCS2)

.PHONY: rosenbrock woods all clean

rosenbrock: $(EXEC1)

woods: $(EXEC2)

all: rosenbrock woods
	$(MV) $(MVFLAGS) $(EXEC1) $(EXEC2) $(BIN_DIR)

clean:
	$(RM) $(EXEC1) $(EXEC2)

# vim:set nu:et:ts=4:sw=4:
