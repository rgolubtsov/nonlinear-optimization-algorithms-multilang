# =============================================================================
# nlp-unconstrained-core/hooke-jeeves/f77/src/Makefile
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

EXEC1   = hooke
EXEC2S  = -woods
EXEC2   = $(EXEC1)$(EXEC2S)
DEPS1P  = rosenbrock
DEPS1   = $(DEPS1P).o
DEPS2P  = woods
DEPS2   = $(DEPS2P).o
BIN_DIR = ../bin
OBJ_DIR = ../obj

# Specify flags and other vars here.
CC      = gfortran
CFLAGS  = -Wall
LDFLAGS = -o

MV      = mv
MVFLAGS = -v

# Making the first target (rosenbrock).
$(DEPS1): CFLAGS += -c -o
$(DEPS1): %.o: %.f
	$(CC) $(CFLAGS) $@ $<

$(EXEC1): $(DEPS1)
	$(CC) $(LDFLAGS) $@ $<

# Making the second target (woods).
$(DEPS2): CFLAGS += -DWoods -c -o
$(DEPS2): %.o: %.f
	$(CC) $(CFLAGS) $@ $<

$(EXEC2): $(DEPS2)
	$(CC) $(LDFLAGS) $@ $<

.PHONY: rosenbrock woods all clean

rosenbrock: $(EXEC1)

woods: $(EXEC2)

all: rosenbrock woods
	$(MV) $(MVFLAGS) $(EXEC1) $(EXEC2) $(BIN_DIR)
	$(MV) $(MVFLAGS) $(DEPS1) $(DEPS2) $(OBJ_DIR)

clean:
	$(RM) $(EXEC1) $(DEPS1) $(EXEC2) $(DEPS2)

# =============================================================================
# vim:set nu:et:ts=4:sw=4:
# =============================================================================
