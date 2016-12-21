/*
 * nlp-unconstrained-cli/nelder-mead/c/src/woods.c
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

    s1 = x[INDEX_1] - pow(x[INDEX_0], SQUARE);
    s2 = ACOEFF     -     x[INDEX_0];
    s3 = x[INDEX_1] -     ACOEFF;

    t1 = x[INDEX_3] - pow(x[INDEX_2], SQUARE);
    t2 = ACOEFF     -     x[INDEX_2];
    t3 = x[INDEX_3] -     ACOEFF;

    t4 = s3 + t3;
    t5 = s3 - t3;

    return (SCOEFF * pow(s1, SQUARE) + pow(s2, SQUARE)
          + TCOEFF * pow(t1, SQUARE) + pow(t2, SQUARE)
          + DCOEFF * pow(t4, SQUARE) + pow(t5, SQUARE) / DCOEFF);
}

/* ========================================================================= */
/* vim:set nu:et:ts=4:sw=4:                                                  */
/* ========================================================================= */
