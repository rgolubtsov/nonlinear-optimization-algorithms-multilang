/*
 * nlp-unconstrained-cli/hooke-jeeves/vala/src/rosenbrock.vala
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * Written by Radislav (Radicchio) Golubtsov, 2016
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

namespace CLIHooke {

// Helper constants.
const uint   INDEX_ZERO             =   0;
const uint   INDEX_ONE              =   1;
const double ONE_HUNDRED_POINT_ZERO = 100.0;
const double ONE_POINT_ZERO         =   1.0;

/**
 * The user-supplied objective function f(x,n).
 * Represents here the Rosenbrock's classic parabolic valley
 * (&quot;banana&quot;) function.
 *
 * @param x        The point at which f(x) should be evaluated.
 * @param n        The number of coordinates of <code>x</code>.
 * @param funevals The number of function evaluations.
 *
 * @return The objective function value.
 */
double f(double *x, uint n, uint funevals) {
    double a;
    double b;
    double c;

    funevals++;

    a = x[INDEX_ZERO];
    b = x[INDEX_ONE];

    c = ONE_HUNDRED_POINT_ZERO * (b - (a * a)) * (b - (a * a));

    return (c + ((ONE_POINT_ZERO - a) * (ONE_POINT_ZERO - a)));
}

} // namespace CLIHooke

// vim:set nu:et:ts=4:sw=4:
