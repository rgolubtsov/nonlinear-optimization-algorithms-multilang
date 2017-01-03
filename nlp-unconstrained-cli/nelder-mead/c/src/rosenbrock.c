/*
 * nlp-unconstrained-cli/nelder-mead/c/src/rosenbrock.c
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Nelder-Mead nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * Written by Radislav (Radicchio) Golubtsov, 2015-2017
 *
 * This is free and unencumbered software released into the public domain.
 *
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 *
 * (See the LICENSE file at the top of the source tree.)
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

/* vim:set nu:et:ts=4:sw=4: */
