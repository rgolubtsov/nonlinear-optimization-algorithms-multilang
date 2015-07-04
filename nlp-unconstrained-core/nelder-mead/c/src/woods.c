/*
 * nlp-unconstrained-core/nelder-mead/c/src/woods.c
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

#include "woods.h"

/* The user-supplied objective function f(x). */
double f(const double *x) {
    double s1;
    double s2;
    double s3;
    double t1;
    double t2;
    double t3;
    double t4;
    double t5;

    s1 = x[1] - pow(x[0], 2);
    s2 = 1    -     x[0];
    s3 = x[1] -     1;

    t1 = x[3] - pow(x[2], 2);
    t2 = 1    -     x[2];
    t3 = x[3] -     1;

    t4 = s3 + t3;
    t5 = s3 - t3;

    return (100 * pow(s1, 2) + pow(s2, 2)
           + 90 * pow(t1, 2) + pow(t2, 2)
           + 10 * pow(t4, 2) + pow(t5, 2) / 10);
}

/* ========================================================================= */
/* vim:set nu:et:ts=4:sw=4:                                                  */
/* ========================================================================= */
