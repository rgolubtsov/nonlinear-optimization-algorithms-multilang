/*
 * nlp-unconstrained-cli/hooke-jeeves/go/src/woods.go
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * Written by Radislav (Radicchio) Golubtsov, 2015-2020
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

/** Helper constants. */
const ONE_HUNDRED uint    = 100
const NINETY      uint    =  90
const TEN         uint    =  10
const TEN_POINT   float64 =  10.

/**
 * The <code>Woods</code> structure is responsible for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 * <br />
 * <br />The objective function in this case
 * is the so-called &quot;Woods&quot; function.
 *
 * @author  Radislav (Radicchio) Golubtsov
 * @version 0.1.1
 * @see     Hooke (hooke.go)
 * @since   hooke-jeeves 0.1
 */
type Woods struct { }

/**
 * The user-supplied objective function F(x,n).
 * <br />
 * <br />Woods &ndash; a la More, Garbow &amp; Hillstrom
 * (TOMS algorithm 566).
 *
 * @param x         The point at which F(x) should be evaluated.
 * @param n         The number of coordinates of <code>x</code>.
 * @param cFunEvals The number of function evaluations container
 *                  (*FunEvals).
 *
 * @return The objective function value.
 */
func (w Woods) F(x []float64, n uint, cFunEvals *FunEvals) float64 {
    var s1 float64
    var s2 float64
    var s3 float64
    var t1 float64
    var t2 float64
    var t3 float64
    var t4 float64
    var t5 float64

    cFunEvals.SetFunEvals(cFunEvals.GetFunEvals() + 1)

    s1 = x[INDEX_ONE]   - x[INDEX_ZERO] * x[INDEX_ZERO]
    s2 = 1              - x[INDEX_ZERO]
    s3 = x[INDEX_ONE]   - 1

    t1 = x[INDEX_THREE] - x[INDEX_TWO]  * x[INDEX_TWO]
    t2 = 1              - x[INDEX_TWO]
    t3 = x[INDEX_THREE] - 1

    t4 = s3 + t3
    t5 = s3 - t3

    return (float64(ONE_HUNDRED) * (s1 * s1) + s2 * s2 +
                 float64(NINETY) * (t1 * t1) + t2 * t2 +
                    float64(TEN) * (t4 * t4) + t5 * t5 / TEN_POINT)
}

// vim:set nu et ts=4 sw=4:
