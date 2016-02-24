/*
 * nlp-unconstrained-core/hooke-jeeves/go/src/rosenbrock.go
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

import "fmt"

/** Constant. The stepsize geometric shrink. */
const RHO_BEGIN float32 = 0.5

/** Helper constants. */
const ONE_HUNDRED_POINT_ZERO float32 = 100.0
const ONE_POINT_ZERO         float32 =   1.0
const TWO                    uint    =   2
const MINUS_ONE_POINT_TWO    float32 =  -1.2

/**
 * The <code>Rosenbrock</code> interface is responsible for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 * <br />
 * <br />The objective function in this case
 * is the Rosenbrock's parabolic valley function.
 *
 * @author  Radislav (Radic) Golubtsov
 * @version 0.1
 * @see     Hooke
 * @since   hooke-jeeves 0.1
 */
type Rosenbrock interface {
    // The user-supplied objective function f(x,n).
    f() float32
}

type sRosenbrock struct { }

/**
 * The user-supplied objective function f(x,n).
 * <br />
 * <br />Represents here the Rosenbrock's classic parabolic valley
 * (&quot;banana&quot;) function.
 *
 * @param x The point at which f(x) should be evaluated.
 * @param n The number of coordinates of <code>x</code>.
 *
 * @return The objective function value.
 */
func (s_rosenbrock sRosenbrock) f(x []float32, n uint) float32 {
    const INDEX_ZERO uint = 0
    const INDEX_ONE  uint = 1

    var a float32
    var b float32
    var c float32

    setFunEvals(getFunEvals() + 1)

    a = x[INDEX_ZERO]
    b = x[INDEX_ONE]

    c = ONE_HUNDRED_POINT_ZERO * (b - (a * a)) * (b - (a * a))

    return (c + ((ONE_POINT_ZERO - a) * (ONE_POINT_ZERO - a)))
}

// Main program function main() :-).
func main() {
    var jj uint

    fmt.Printf("\n\n\nHOOKE USED %d ITERATIONS, AND RETURNED\n", jj)
}

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
