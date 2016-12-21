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

    a = ACOEFF     -     x[INDEX_0];
    b = x[INDEX_1] - pow(x[INDEX_0], SQUARE);

    return (pow(a, SQUARE) + BCOEFF * pow(b, SQUARE));
}

/* ========================================================================= */
/* vim:set nu:et:ts=4:sw=4:                                                  */
/* ========================================================================= */
