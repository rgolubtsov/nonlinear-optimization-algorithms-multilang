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

import "reflect"

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

/**
 * Helper method.
 * <br />
 * <br />Given a point, look for a better one nearby, one coord at a time.
 *
 * @param delta     The delta between <code>prevBest</code>
 *                  and <code>point</code>.
 * @param point     The coordinate from where to begin.
 * @param prevBest  The previous best-valued coordinate.
 * @param nVars     The number of variables.
 * @param objFunCls The class in which the objective function is defined.
 *
 * @return The objective function value at a nearby.
 */
func bestNearby(delta []float32, point []float32, prevBest float32, nVars uint, objFunCls struct{}) float32 {
        var minF       float32
        var z    [VARS]float32
        var fTmp       float32

        var i uint

        var s_rosenbrock sRosenbrock
        var s_woods      sWoods

        minF = prevBest

        for i = 0; i < nVars; i++ {
            z[i] = point[i]
        }

        for i = 0; i < nVars; i++ {
            z[i] = point[i] + delta[i]

            if reflect.DeepEqual(objFunCls, s_rosenbrock) {
//                fTmp = (*sRosenbrock)(nil).f(z[0:], nVars)
                fTmp = s_rosenbrock.f(z[0:], nVars)
            } else if reflect.DeepEqual(objFunCls, s_woods) {
//                fTmp = (*sWoods)(nil).f(z[0:], nVars)
                fTmp = s_woods.f(z[0:], nVars)
            } else {
                fTmp = 0
            }

            if fTmp < minF {
                minF = fTmp
            } else {
                delta[i] = 0.0 - delta[i]
                z[i]     = point[i] + delta[i]

                if reflect.DeepEqual(objFunCls, s_rosenbrock) {
//                    fTmp = (*sRosenbrock)(nil).f(z[0:], nVars)
                    fTmp = s_rosenbrock.f(z[0:], nVars)
                } else if reflect.DeepEqual(objFunCls, s_woods) {
//                    fTmp = (*sWoods)(nil).f(z[0:], nVars)
                    fTmp = s_woods.f(z[0:], nVars)
                } else {
                    fTmp = 0
                }

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

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
