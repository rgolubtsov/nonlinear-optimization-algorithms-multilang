/*
 * nlp-unconstrained-core/hooke-jeeves/go/src/funevals.go
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

/**
 * The <code>FunEvals</code> structure is a helper structure.
 * It holds the only property &ndash; the number of objective function
 * evaluations (and corresponding accessor methods).
 *
 * @author  Radislav (Radic) Golubtsov
 * @version 0.1
 * @see     Hooke (hooke.go)
 * @since   hooke-jeeves 0.1
 */
type FunEvals struct {
    /** The number of function evaluations. */
    funEvals uint
}

/**
 * Getter for <code>funEvals</code>.
 *
 * @return The number of function evaluations.
 */
func (fe *FunEvals) GetFunEvals() uint {
    return fe.funEvals
}

/**
 * Setter for <code>funEvals</code>.
 *
 * @param __funEvals The number of function evaluations.
 */
func (fe *FunEvals) SetFunEvals(__funEvals uint) {
    fe.funEvals = __funEvals
}

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
