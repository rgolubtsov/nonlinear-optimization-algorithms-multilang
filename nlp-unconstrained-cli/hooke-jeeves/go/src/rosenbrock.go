/*
 * nlp-unconstrained-cli/hooke-jeeves/go/src/rosenbrock.go
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * Written by Radislav (Radicchio) Golubtsov, 2015-2025
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

package main

/** Helper constant. */
const ONE_HUNDRED_POINT_ZERO float64 = 100.0

/**
 * The <code>Rosenbrock</code> structure is responsible for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 * <br />
 * <br />The objective function in this case
 * is the Rosenbrock's parabolic valley function.
 *
 * @author  Radislav (Radicchio) Golubtsov
 * @version 0.1.1
 * @see     Hooke (hooke.go)
 * @since   hooke-jeeves 0.1
 */
type Rosenbrock struct { }

/**
 * The user-supplied objective function F(x,n).
 * <br />
 * <br />Represents here the Rosenbrock's classic parabolic valley
 * (&quot;banana&quot;) function.
 *
 * @param x         The point at which F(x) should be evaluated.
 * @param n         The number of coordinates of <code>x</code>.
 * @param cFunEvals The number of function evaluations container
 *                  (*FunEvals).
 *
 * @return The objective function value.
 */
func (r Rosenbrock) F(x []float64, n uint, cFunEvals *FunEvals) float64 {
    var a float64
    var b float64
    var c float64

    cFunEvals.SetFunEvals(cFunEvals.GetFunEvals() + 1)

    a = x[INDEX_ZERO]
    b = x[INDEX_ONE]

    c = ONE_HUNDRED_POINT_ZERO * (b - (a * a)) * (b - (a * a))

    return (c + ((ONE_POINT_ZERO - a) * (ONE_POINT_ZERO - a)))
}

// vim:set nu et ts=4 sw=4:
