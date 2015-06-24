# =============================================================================
# nlp-unconstrained-core/hooke-jeeves/__orig/src/Makefile
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
# =============================================================================

EXEC1   = hooke
EXEC2S  = -woods
EXEC2   = $(EXEC1)$(EXEC2S)
DEPS1   = $(EXEC1).o
DEPS2   = $(EXEC2).o
BIN_DIR = ../bin
OBJ_DIR = ../obj

# Specify flags and other vars here.
# Note: To use the system default C compiler (likely gcc, the GNU C Compiler)
#       instead of tcc (Tiny C Compiler), simply comment out the next line,
#       just like this: "#CC = tcc" - and it will be used automatically.
#CC     = tcc
CFLAGS = -Wall -pedantic -O3 -Wtraditional -Wno-implicit-int -Wno-return-type

# -lm - (Math library) - Needs to be specified explicitly
#                        only when linking through tcc.
#LDLIBS = -lm

MV      = mv
MVFLAGS = -v

# Making the first target (rosenbrock).
$(DEPS1): %.o: %.c

$(EXEC1): $(DEPS1)

# Making the second target (woods).
$(DEPS2): CFLAGS += -DWOODS -c -o
$(DEPS2): %$(EXEC2S).o: %.c
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

# =============================================================================
# vim:set nu:et:ts=4:sw=4:
# =============================================================================
