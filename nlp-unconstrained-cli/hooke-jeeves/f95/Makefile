#
# nlp-unconstrained-cli/hooke-jeeves/f95/Makefile
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
# =============================================================================
# Written by Radislav (Radicchio) Golubtsov, 2015-2024
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

BIN_DIR    = bin
OBJ_DIR    = obj
SRC_DIR    = src
ALL_TARGET = all

# Specify flags and other vars here.
MKDIR      = mkdir

# Using here this user-defined var because "The options ‘-C’, ‘-f’, ‘-o’,
# and ‘-W’ are not put into MAKEFLAGS; these options are not passed down."
# -- from the GNU Make Manual. I.e. use of MAKEFLAGS cannot be applicable
# in this case.
MAKE_FLAGS = -C

RMFLAGS    = -vR

# Making the target.
$(BIN_DIR) $(OBJ_DIR):
	if [ ! -d "$(BIN_DIR)" ]; then \
		$(MKDIR) $(BIN_DIR); \
		if [ ! -d "$(OBJ_DIR)" ]; then \
			$(MKDIR) $(OBJ_DIR); \
		fi; \
		$(MAKE) $(MAKE_FLAGS)$(SRC_DIR) $(ALL_TARGET); \
	fi

.PHONY: all clean

all: $(BIN_DIR) $(OBJ_DIR)

clean:
	$(RM) $(RMFLAGS) $(BIN_DIR) $(OBJ_DIR)

# vim:set nu ts=4 sw=4:
