/*
 * nlp-unconstrained-core/hooke-jeeves/go/src/woods.go
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

/**
 * Constant. The stepsize geometric shrink.
 * <br />
 * <br />The Hooke &amp; Jeeves algorithm works reasonably well
 * on Rosenbrock's function, but can fare worse on some standard
 * test functions, depending on rho. Here is an example that works well
 * when rho = 0.5, but fares poorly with rho = 0.6, and better again
 * with rho = 0.8.
 */
const RHO_WOODS float32 = 0.6

/** Helper constants. */
const INDEX_TWO   uint    =  2
const INDEX_THREE uint    =  3
const ONE_HUNDRED uint    =  100
const NINETY      uint    =  90
const TEN         uint    =  10
const TEN_POINT   float32 =  10.
const FOUR        uint    =  4
const MINUS_THREE int     = -3
const MINUS_ONE   int     = -1

/**
 * The <code>Woods</code> interface is responsible for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 * <br />
 * <br />The objective function in this case
 * is the so-called &quot;Woods&quot; function.
 *
 * @author  Radislav (Radic) Golubtsov
 * @version 0.1
 * @see     Hooke
 * @since   hooke-jeeves 0.1
 */
type Woods interface {
    // The user-supplied objective function f(x,n).
    f() float32
}

type sWoods struct { }

/**
 * The user-supplied objective function f(x,n).
 * <br />
 * <br />Woods &ndash; a la More, Garbow &amp; Hillstrom
 * (TOMS algorithm 566).
 *
 * @param x The point at which f(x) should be evaluated.
 * @param n The number of coordinates of <code>x</code>.
 *
 * @return The objective function value.
 */
func (s_woods sWoods) f(x []float32, n uint) float32 {
    const INDEX_ONE  uint = 1
    const INDEX_ZERO uint = 0

    var s1 float32
    var s2 float32
    var s3 float32
    var t1 float32
    var t2 float32
    var t3 float32
    var t4 float32
    var t5 float32

    setFunEvals(getFunEvals() + 1)

    s1 = x[INDEX_ONE]   - x[INDEX_ZERO] * x[INDEX_ZERO]
    s2 = 1              - x[INDEX_ZERO]
    s3 = x[INDEX_ONE]   - 1

    t1 = x[INDEX_THREE] - x[INDEX_TWO]  * x[INDEX_TWO]
    t2 = 1              - x[INDEX_TWO]
    t3 = x[INDEX_THREE] - 1

    t4 = s3 + t3
    t5 = s3 - t3

    return (float32(ONE_HUNDRED) * (s1 * s1) + s2 * s2 +
                 float32(NINETY) * (t1 * t1) + t2 * t2 +
                    float32(TEN) * (t4 * t4) + t5 * t5 / TEN_POINT)
}

// Main program function main() :-).
func main() {
    var jj uint

    fmt.Printf("\n\n\nHOOKE USED %d ITERATIONS, AND RETURNED\n", jj)

    fmt.Println("True answer: f(1, 1, 1, 1) = 0.")
}

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
