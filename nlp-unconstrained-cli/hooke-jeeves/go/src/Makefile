#
# nlp-unconstrained-cli/hooke-jeeves/go/src/Makefile
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
# Written by Radislav (Radicchio) Golubtsov, 2015-2020
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

EXEC    = hooke
DEPS    = $(EXEC).o
SRCS    = $(EXEC).go    \
          rosenbrock.go \
          woods.go      \
          funevals.go
BIN_DIR = ../bin
OBJ_DIR = ../obj

# Specify flags and other vars here.
GO            = go tool
GO_COMP_FLAGS = compile -complete -o
GO_LINK_FLAGS = link -s -w -o

MV      = mv
MVFLAGS = -v

# Making the target.
$(DEPS): $(SRCS)
	$(GO) $(GO_COMP_FLAGS) $(DEPS) $(SRCS)

$(EXEC): $(DEPS)
	$(GO) $(GO_LINK_FLAGS) $(EXEC) $(DEPS)

.PHONY: all clean

all: $(EXEC)
	$(MV) $(MVFLAGS) $(EXEC) $(BIN_DIR)
	$(MV) $(MVFLAGS) $(DEPS) $(OBJ_DIR)

clean:
	$(RM) $(EXEC) $(DEPS)

# vim:set nu ts=4 sw=4:
