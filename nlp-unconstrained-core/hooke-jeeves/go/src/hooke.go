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

import "math"
import "fmt"
import "os"

/** Constant. The maximum number of variables. */
const VARS uint = 250

/** Constant. The stepsize geometric shrink. */
const RHO_BEGIN float64 = 0.5

/**
 * Constant. The stepsize geometric shrink.
 * <br />
 * <br />The Hooke &amp; Jeeves algorithm works reasonably well
 * on Rosenbrock's function, but can fare worse on some standard
 * test functions, depending on rho. Here is an example that works well
 * when rho = 0.5, but fares poorly with rho = 0.6, and better again
 * with rho = 0.8.
 */
const RHO_WOODS float64 = 0.6

/** Constant. The ending value of stepsize. */
const EPSMIN float64 = 1E-6

/** Constant. The maximum number of iterations. */
const IMAX uint = 5000

/** Helper constants. */
const WOODS               string  = "woods"
const INDEX_ZERO          uint    =  0
const INDEX_ONE           uint    =  1
const INDEX_TWO           uint    =  2
const INDEX_THREE         uint    =  3
const TWO                 uint    =  2
const FOUR                uint    =  4
const MINUS_ONE_POINT_TWO float64 = -1.2
const ONE_POINT_ZERO      float64 =  1.0
const MINUS_THREE         int     = -3
const MINUS_ONE           int     = -1
const ZERO_POINT_FIVE     float64 =  0.5

/**
 * The <code>Hooke</code> structure contains methods for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 *
 * @author  Radislav (Radic) Golubtsov
 * @version 0.1
 * @see     Rosenbrock (rosenbrock.go)
 * @see     Woods           (woods.go)
 * @since   hooke-jeeves 0.1
 */
type Hooke struct {
    /** The number of function evaluations. */
    funEvals uint
}

/**
 * Getter for <code>funEvals</code>.
 *
 * @return The number of function evaluations.
 */
func (h Hooke) GetFunEvals() uint {
    return h.funEvals
}

/**
 * Setter for <code>funEvals</code>.
 *
 * @param __funEvals The number of function evaluations.
 */
func (h Hooke) SetFunEvals(__funEvals uint) {
    h.funEvals = __funEvals
}

/**
 * Helper method.
 * <br />
 * <br />Given a point, look for a better one nearby, one coord at a time.
 *
 * @param delta    The delta between <code>prevBest</code>
 *                 and <code>point</code>.
 * @param point    The coordinate from where to begin.
 * @param prevBest The previous best-valued coordinate.
 * @param nVars    The number of variables.
 * @param woods    The cmd-line arg to indicate,
 *                 which objective function to use.
 *
 * @return The objective function value at a nearby.
 */
func (h Hooke) BestNearby(delta    []float64,
                          point    []float64,
                          prevBest   float64,
                          nVars      uint,
                          woods      string) float64 {

    var minF       float64
    var z    [VARS]float64
    var fTmp       float64

    var i uint

    minF = prevBest

    for i = 0; i < nVars; i++ {
        z[i] = point[i]
    }

    r := new(Rosenbrock)
    w := new(Woods)

    for i = 0; i < nVars; i++ {
        z[i] = point[i] + delta[i]

        if woods != WOODS {                                    // #ifndef WOODS
            fTmp = r.F(z[0:], nVars)
        } else {                                               // #else
            fTmp = w.F(z[0:], nVars)
        }                                                      // #endif

        if fTmp < minF {
            minF = fTmp
        } else {
            delta[i] = 0.0 - delta[i]
            z[i]     = point[i] + delta[i]

            if woods != WOODS {                                // #ifndef WOODS
                fTmp = r.F(z[0:], nVars)
            } else {                                           // #else
                fTmp = w.F(z[0:], nVars)
            }                                                  // #endif

            if fTmp < minF {
                minF = fTmp
            } else {
                z[i] = point[i]
            }
        }
    }

    for i = 0; i < nVars; i++ {
        point[i] = z[i]
    }

    return minF
}

/**
 * Main optimization method.
 * <br />
 * <br />The hooke subroutine itself.
 *
 * @param nVars   The number of variables.
 * @param startPt The starting point coordinates.
 * @param endPt   The ending point coordinates.
 * @param rho     The rho value.
 * @param epsilon The epsilon value.
 * @param iterMax The maximum number of iterations.
 * @param woods   The cmd-line arg to indicate,
 *                which objective function to use.
 *
 * @return The number of iterations used to find the local minimum.
 */
func (h Hooke) hooke(nVars     uint,
                     startPt []float64,
                     endPt   []float64,
                     rho       float64,
                     epsilon   float64,
                     iterMax   uint,
                     woods     string) uint {

    var i     uint
    var iAdj  uint
    var iters uint
    var j     uint
    var keep  uint

    var newX       [VARS]float64
    var xBefore    [VARS]float64
    var delta      [VARS]float64
    var stepLength       float64
    var fBefore          float64
    var newF             float64
    var tmp              float64

    for i = 0; i < nVars; i++ {
        xBefore[i] = startPt[i]
        newX[i]    = xBefore[i]

        delta[i] = math.Abs(startPt[i] * rho)

        if delta[i] == 0.0 {
            delta[i] = rho
        }
    }

    iAdj       = 0
    stepLength = rho
    iters      = 0

    r := new(Rosenbrock)
    w := new(Woods)

    if woods != WOODS {                                        // #ifndef WOODS
        fBefore = r.F(newX[0:], nVars)
    } else {                                                   // #else
        fBefore = w.F(newX[0:], nVars)
    }                                                          // #endif

    newF = fBefore

    for (iters < iterMax) && (stepLength > epsilon) {
        iters++
        iAdj++

        fmt.Printf(
            "\nAfter %5d funevals, f(x) =  %.4e at\n", h.funEvals, fBefore)

        for j = 0; j < nVars; j++ {
            fmt.Printf("   x[%2d] = %.4e\n", j, xBefore[j])
        }

        // Find best new point, one coord at a time.
        for i = 0; i < nVars; i++ {
            newX[i] = xBefore[i]
        }

        newF = h.BestNearby(delta[0:], newX[0:], fBefore, nVars, woods)

        // If we made some improvements, pursue that direction.
        keep = 1

        for (newF < fBefore) && (keep == 1) {
            iAdj = 0

            for i = 0; i < nVars; i++ {
                // Firstly, arrange the sign of delta[].
                if newX[i] <= xBefore[i] {
                    delta[i] = 0.0 - math.Abs(delta[i])
                } else {
                    delta[i] = math.Abs(delta[i])
                }

                // Now, move further in this direction.
                tmp        = xBefore[i]
                xBefore[i] = newX[i]
                newX[i]    = newX[i] + newX[i] - tmp
            }

            fBefore = newF

            newF = h.BestNearby(delta[0:], newX[0:], fBefore, nVars, woods)

            // If the further (optimistic) move was bad....
            if newF >= fBefore {
                break
            }

            /*
             * Make sure that the differences between the new and the old
             * points are due to actual displacements; beware of roundoff
             * errors that might cause newF < fBefore.
             */
            keep = 0

            for i = 0; i < nVars; i++ {
                keep = 1

                if math.Abs(newX[i] - xBefore[i]) >
                    (ZERO_POINT_FIVE * math.Abs(delta[i])) {

                    break
                } else {
                    keep = 0
                }
            }
        }

        if (stepLength >= epsilon) && (newF >= fBefore) {
            stepLength = stepLength * rho

            for i = 0; i < nVars; i++ {
                delta[i] *= rho
            }
        }
    }

    for i = 0; i < nVars; i++ {
        endPt[i] = xBefore[i]
    }

    return iters
}

/*
 * Main program function main() :-).
 *
 * It looks for the presence of the 'woods' cmd-line arg as the first arg,
 * and might be run without any args -- Rosenbrock test problem will be solved:
 *
 *     $ ./nlp-unconstrained-core/hooke-jeeves/go/bin/hooke
 *
 * Or it might be run with the 'woods' arg --Woods test problem will be solved:
 *
 *     $ ./nlp-unconstrained-core/hooke-jeeves/go/bin/hooke woods
 */
func main() {
    var nVars   uint
    var iterMax uint
    var jj      uint
    var i       uint

    var startPt [VARS]float64
    var rho           float64
    var epsilon       float64
    var endPt   [VARS]float64

    var argsLen uint = uint(len(os.Args) - 1)

    var woods string

    if argsLen > 0 {
        woods = os.Args[1]
    }

    if woods != WOODS {                                        // #ifndef WOODS
        // Starting guess for Rosenbrock's test function.
        nVars                = TWO
        startPt[INDEX_ZERO]  = MINUS_ONE_POINT_TWO
        startPt[INDEX_ONE]   = ONE_POINT_ZERO
        rho                  = RHO_BEGIN
    } else {                                                   // #else
        // Starting guess test problem "Woods".
        nVars                = FOUR
        startPt[INDEX_ZERO]  = float64(MINUS_THREE)
        startPt[INDEX_ONE]   = float64(MINUS_ONE)
        startPt[INDEX_TWO]   = float64(MINUS_THREE)
        startPt[INDEX_THREE] = float64(MINUS_ONE)
        rho                  = RHO_WOODS
    }                                                          // #endif

    iterMax = IMAX
    epsilon = EPSMIN

    h := new(Hooke)

    jj = h.hooke(nVars, startPt[0:], endPt[0:], rho, epsilon, iterMax, woods)

    fmt.Printf("\n\n\nHOOKE USED %d ITERATIONS, AND RETURNED\n", jj)

    for i = 0; i < nVars; i++ {
        fmt.Printf("x[%3d] = %15.7e \n", i, endPt[i])
    }

    if woods == WOODS {                                         // #ifdef WOODS
        fmt.Println("True answer: f(1, 1, 1, 1) = 0.")
    }                                                           // #endif
}

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
