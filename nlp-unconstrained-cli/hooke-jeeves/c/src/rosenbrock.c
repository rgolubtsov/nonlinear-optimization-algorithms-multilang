/*
 * nlp-unconstrained-cli/hooke-jeeves/c/src/rosenbrock.c
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

#include "rosenbrock.h"
#include "funevals.h"

/* The user-supplied objective function f(x,n). */
double f(const double *x, const unsigned int n, void *__fun_evals) {
    double a;
    double b;
    double c;

    set_funevals(__fun_evals, get_funevals(__fun_evals) + 1);

    a = x[INDEX_ZERO];
    b = x[INDEX_ONE];

    c = ONE_HUNDRED_POINT_ZERO * (b - (a * a)) * (b - (a * a));

    return (c + ((ONE_POINT_ZERO - a) * (ONE_POINT_ZERO - a)));
}

/* ========================================================================= */
/* vim:set nu:et:ts=4:sw=4:                                                  */
/* ========================================================================= */
