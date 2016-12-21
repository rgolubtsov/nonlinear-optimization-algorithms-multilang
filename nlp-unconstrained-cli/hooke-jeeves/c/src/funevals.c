/*
 * nlp-unconstrained-cli/hooke-jeeves/c/src/funevals.c
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 */

#include "funevals.h"

/* Getter for fun_evals.funevals. */
unsigned int get_funevals(const struct fun_evals *__fun_evals) {
    return __fun_evals->funevals;
}

/* Setter for fun_evals.funevals. */
void set_funevals(struct fun_evals *__fun_evals,
                  const unsigned int __funevals) {

    __fun_evals->funevals = __funevals;
}

/* ========================================================================= */
/* vim:set nu:et:ts=4:sw=4:                                                  */
/* ========================================================================= */
