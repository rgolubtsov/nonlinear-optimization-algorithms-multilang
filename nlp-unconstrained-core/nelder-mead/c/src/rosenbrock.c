/*
 * nlp-unconstrained-core/nelder-mead/c/src/rosenbrock.c
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Nelder-Mead nonlinear unconstrained minimization algorithm.
 * ============================================================================
 */

#include "rosenbrock.h"

/* The user-supplied objective function f(x). */
double f(const double *x) {
    double a;
    double b;

    a = x[1] - pow(x[0], 2);
    b = 1    -     x[0];

    return (100 * pow(a, 2) + pow(b, 2));
}

/* ========================================================================= */
/* vim:set nu:et:ts=4:sw=4:                                                  */
/* ========================================================================= */
