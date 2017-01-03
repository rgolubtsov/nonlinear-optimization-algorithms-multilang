#
# nlp-unconstrained-cli/hooke-jeeves/perl/Makefile
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
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

BIN_DIR      = bin
LIB_DIR      = lib
SRC_DIR      = src
BASE_PACKAGE = NLPUCCLIHooke
BASE_CLASS   = $(BASE_PACKAGE).pm
CLASSES      = $(SRC_DIR)/$(BASE_CLASS) $(SRC_DIR)/$(BASE_PACKAGE)
SCRIPT1      = hooke
SCRIPT2S     = -woods
SCRIPT2      = $(SCRIPT1)$(SCRIPT2S)

# Specify flags and other vars here.
PERL      = perl
PERLFLAGS = -c

CD           = cd
GREP         = grep
EXIT         = exit
EXIT_FAILURE = 1
EXIT_SUCCESS = 0
MKDIR        = mkdir
ECHO         = echo
ECHOFLAGS    = -e
SED          = sed
SEDFLAGS     = -i -e
SEDREGEX     = "s/^\ //g"
CHMOD        = chmod
CHMODFLAGS   = -v
CHMODMODE    = 755
CP           = cp
CPFLAGS      = -vR
RMFLAGS      = -vR

# Startup scripts for Rosenbrock and Woods test problems like here documents.
# Rosenbrock:
define SCRIPT1_FULL
#!/usr/bin/env perl\n\
# =============================================================================\n\
# nlp-unconstrained-cli/hooke-jeeves/perl/$(BIN_DIR)/$(SCRIPT1)\n\
# =============================================================================\n\
# Nonlinear Optimization Algorithms Multilang. Version 0.1\n\
# =============================================================================\n\
\n\
use strict;\n\
use warnings;\n\
use v5.10;\n\
\n\
use lib \"$(LIB_DIR)\";\n\
\n\
use $(BASE_PACKAGE);\n\
\n\
# Solving the Rosenbrock test problem.\n\
$(BASE_PACKAGE)->rosenbrock();\n\
\n\
# =============================================================================\n\
# vim:set nu:et:ts=4:sw=4:\n\
# =============================================================================
endef

# Woods:
define SCRIPT2_FULL
#!/usr/bin/env perl\n\
# =============================================================================\n\
# nlp-unconstrained-cli/hooke-jeeves/perl/$(BIN_DIR)/$(SCRIPT2)\n\
# =============================================================================\n\
# Nonlinear Optimization Algorithms Multilang. Version 0.1\n\
# =============================================================================\n\
\n\
use strict;\n\
use warnings;\n\
use v5.10;\n\
\n\
use lib \"$(LIB_DIR)\";\n\
\n\
use $(BASE_PACKAGE);\n\
\n\
# Solving the Woods test problem.\n\
$(BASE_PACKAGE)->woods();\n\
\n\
# =============================================================================\n\
# vim:set nu:et:ts=4:sw=4:\n\
# =============================================================================
endef

# Making the target.
$(BIN_DIR) $(LIB_DIR):
	$(CD) $(SRC_DIR); \
	if [ -n "$$($(PERL) $(PERLFLAGS) $(BASE_CLASS) 2>&1 | $(GREP) 'failed')" ]; then \
		$(ECHO) "Error: Something gets wrong... Exiting..."; \
		$(EXIT) $(EXIT_FAILURE); \
	fi; \
	$(CD) -

	if [ ! -d "$(BIN_DIR)" ]; then \
		$(MKDIR) $(BIN_DIR); \
		if [ ! -d "$(LIB_DIR)" ]; then \
			$(MKDIR) $(LIB_DIR); \
		fi; \
		$(ECHO) $(ECHOFLAGS) "$(SCRIPT1_FULL)" > $(BIN_DIR)/$(SCRIPT1); \
		$(ECHO) $(ECHOFLAGS) "$(SCRIPT2_FULL)" > $(BIN_DIR)/$(SCRIPT2); \
		$(SED) $(SEDFLAGS) $(SEDREGEX) $(BIN_DIR)/*; \
		$(CHMOD) $(CHMODFLAGS) $(CHMODMODE) $(BIN_DIR)/*; \
		$(CP) $(CPFLAGS) $(CLASSES) $(LIB_DIR); \
	fi

.PHONY: all clean

all: $(BIN_DIR) $(LIB_DIR)

clean:
	$(RM) $(RMFLAGS) $(BIN_DIR) $(LIB_DIR)

# vim:set nu:et:ts=4:sw=4:
