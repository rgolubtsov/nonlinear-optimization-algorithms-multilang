# =============================================================================
# nlp-unconstrained-core/hooke-jeeves/c/src/Makefile
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
DEPS1S  = rosenbrock
DEPS1   = $(EXEC1).o $(DEPS1S).o
DEPS2S  = woods
DEPS2   = $(EXEC2).o $(DEPS2S)$(EXEC2S).o
DEPSN1  = funevals
DEPSN   = $(DEPSN1).o
BIN_DIR = ../bin
OBJ_DIR = ../obj

# Specify flags and other vars here.
# Note: To use the system default C compiler (likely gcc, the GNU C Compiler)
#       instead of tcc (Tiny C Compiler), simply comment out the next line,
#       just like this: "#CC = tcc" - and it will be used automatically.
#CC     = tcc

#C_STD = c90 # Same as c89, same as -ansi .
#C_STD = c99
C_STD = c11

CFLAGS = -Wall -pedantic -O2 -std=$(C_STD) # The -std=XXX option is not allowed
                                           # when using tcc.
#LDLIBS = -lm                              # <== Needs only when linking
                                           #     through tcc.

MV      = mv
MVFLAGS = -v

# Making the first target (rosenbrock).
$(DEPS1) $(DEPSN): %.o: %.c

$(EXEC1): $(DEPS1) $(DEPSN)

# Making the second target (woods).
$(DEPS2): CFLAGS += -DWOODS -c -o
$(DEPS2): %$(EXEC2S).o: %.c
	$(CC) $(CFLAGS) $@ $<

$(EXEC2): $(DEPS2) $(DEPSN)

.PHONY: rosenbrock woods all clean

rosenbrock: $(EXEC1)

woods: $(EXEC2)

all: rosenbrock woods
	$(MV) $(MVFLAGS) $(EXEC1) $(EXEC2) $(BIN_DIR)
	$(MV) $(MVFLAGS) $(DEPS1) $(DEPS2) $(DEPSN) $(OBJ_DIR)

clean:
	$(RM) $(EXEC1) $(DEPS1) $(EXEC2) $(DEPS2) $(DEPSN)

# =============================================================================
# vim:set nu:et:ts=4:sw=4:
# =============================================================================
