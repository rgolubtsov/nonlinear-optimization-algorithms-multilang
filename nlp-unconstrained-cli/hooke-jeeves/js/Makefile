#
# nlp-unconstrained-cli/hooke-jeeves/js/Makefile
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
# =============================================================================
# Written by Radislav (Radicchio) Golubtsov, 2015-2019
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

BIN_DIR  = bin
LIB_DIR  = lib
SRC_DIR  = src
CLASS1   = rosenbrock.js
CLASS2   = woods.js
CLASS3   = hooke.js
CLASSES  = $(SRC_DIR)/$(CLASS1) $(SRC_DIR)/$(CLASS2) $(SRC_DIR)/$(CLASS3)
SCRIPT1  = hooke
SCRIPT2S = -woods
SCRIPT2  = $(SCRIPT1)$(SCRIPT2S)

# Specify flags and other vars here.
# Node.js:
NODEJS = node

MKDIR      = mkdir
ECHO       = echo
ECHOFLAGS  = -e
SED        = sed
SEDFLAGS   = -i -e
SEDREGEX   = "s/^\ //g;s/'//g"
CHMOD      = chmod
CHMODFLAGS = -v
CHMODMODE  = 755
CP         = cp
CPFLAGS    = -v
RMFLAGS    = -vR

# Startup scripts for Rosenbrock and Woods test problems like here documents.
# Rosenbrock:
define SCRIPT1_FULL
#!/usr/bin/env sh\n\
# =============================================================================\n\
# nlp-unconstrained-cli/hooke-jeeves/js/$(BIN_DIR)/$(SCRIPT1)\n\
# =============================================================================\n\
# Nonlinear Optimization Algorithms Multilang. Version 0.1\n\
# =============================================================================\n\
\n\
declare -r EXIT_FAILURE=1 # Failing exit status.\n\
declare -r EXIT_SUCCESS=0 # Successful exit status.\n\
\n\
declare -r NODEJS=$(NODEJS)\n\
declare -r LIB_DIR=$(LIB_DIR)\n\
declare -r CLASS=$(CLASS1)\n\
\n\
$$'{NODEJS}' $$'{LIB_DIR}'/$$'{CLASS}'\n\
\n\
if [ $$'?' -eq $$'{EXIT_FAILURE}' ]; then\n\
    exit $$'{EXIT_FAILURE}'\n\
fi\n\
\n\
exit $$'{EXIT_SUCCESS}'\n\
\n\
# =============================================================================\n\
# vim:set nu et ts=4 sw=4:\n\
# =============================================================================
endef

# Woods:
define SCRIPT2_FULL
#!/usr/bin/env sh\n\
# =============================================================================\n\
# nlp-unconstrained-cli/hooke-jeeves/js/$(BIN_DIR)/$(SCRIPT2)\n\
# =============================================================================\n\
# Nonlinear Optimization Algorithms Multilang. Version 0.1\n\
# =============================================================================\n\
\n\
declare -r EXIT_FAILURE=1 # Failing exit status.\n\
declare -r EXIT_SUCCESS=0 # Successful exit status.\n\
\n\
declare -r NODEJS=$(NODEJS)\n\
declare -r LIB_DIR=$(LIB_DIR)\n\
declare -r CLASS=$(CLASS2)\n\
\n\
$$'{NODEJS}' $$'{LIB_DIR}'/$$'{CLASS}'\n\
\n\
if [ $$'?' -eq $$'{EXIT_FAILURE}' ]; then\n\
    exit $$'{EXIT_FAILURE}'\n\
fi\n\
\n\
exit $$'{EXIT_SUCCESS}'\n\
\n\
# =============================================================================\n\
# vim:set nu et ts=4 sw=4:\n\
# =============================================================================
endef

# Making the target.
$(BIN_DIR) $(LIB_DIR):
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

# vim:set nu ts=4 sw=4:
