# =============================================================================
# nlp-unconstrained-core/hooke-jeeves/__orig/Makefile
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

# =============================================================================
# vim:set nu:et:ts=4:sw=4:
# =============================================================================
