/*
 * nlp-unconstrained-core/hooke-jeeves/go/src/hooke.go
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

package main

/** Constant. The maximum number of variables. */
const VARS uint = 250

/** Constant. The ending value of stepsize. */
const EPSMIN float32 = 1E-6

/** Constant. The maximum number of iterations. */
const IMAX uint = 5000

/** Helper constants. */
const INDEX_ZERO      uint    = 0
const INDEX_ONE       uint    = 1
const ZERO_POINT_FIVE float32 = 0.5

/** The number of function evaluations. */
var funEvals uint = 0

/**
 * Getter for <code>funEvals</code>.
 *
 * @return The number of function evaluations.
 */
func getFunEvals() uint {
    return funEvals
}

/**
 * Setter for <code>funEvals</code>.
 *
 * @param __funEvals The number of function evaluations.
 */
func setFunEvals(__funEvals uint) {
    funEvals = __funEvals
}

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
