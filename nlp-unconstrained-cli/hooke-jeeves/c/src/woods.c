/*
 * nlp-unconstrained-core/hooke-jeeves/c/src/woods.c
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

#include "woods.h"
#include "funevals.h"

/* The user-supplied objective function f(x,n). */
double f(const double *x, const unsigned int n, void *__fun_evals) {
    double s1;
    double s2;
    double s3;
    double t1;
    double t2;
    double t3;
    double t4;
    double t5;

    set_funevals(__fun_evals, get_funevals(__fun_evals) + 1);

    s1 = x[INDEX_ONE] - x[INDEX_ZERO] * x[INDEX_ZERO];
    s2 = 1 - x[INDEX_ZERO];
    s3 = x[INDEX_ONE] - 1;
    t1 = x[INDEX_THREE] - x[INDEX_TWO] * x[INDEX_TWO];
    t2 = 1 - x[INDEX_TWO];
    t3 = x[INDEX_THREE] - 1;
    t4 = s3 + t3;
    t5 = s3 - t3;

    return (ONE_HUNDRED * (s1 * s1) + s2 * s2
               + NINETY * (t1 * t1) + t2 * t2
                  + TEN * (t4 * t4) + t5 * t5 / TEN_POINT);
}

/* ========================================================================= */
/* vim:set nu:et:ts=4:sw=4:                                                  */
/* ========================================================================= */
