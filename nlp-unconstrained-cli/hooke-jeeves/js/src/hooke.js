/*
 * nlp-unconstrained-cli/hooke-jeeves/js/src/hooke.js
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * Written by Radislav (Radicchio) Golubtsov, 2015-2019
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

"use strict";

/**
 * The <code>NLPUCCLIHooke</code> namespace is used as a container
 * for the <code>Hooke</code> class.
 *
 * @namespace NLPUCCLIHooke
 */
var NLPUCCLIHooke = NLPUCCLIHooke || {};

/**
 * The <code>Hooke</code> class contains methods for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 *
 * @class    Hooke
 * @memberof NLPUCCLIHooke
 * @author   Radislav (Radicchio) Golubtsov
 * @version  0.1
 * @see      {@link https://raw.githubusercontent.com/rgolubtsov/nonlinear-optimization-algorithms-multilang/master/nlp-unconstrained-cli/hooke-jeeves/js/src/rosenbrock.js|NLPUCCLIHooke.Rosenbrock}
 * @see      {@link https://raw.githubusercontent.com/rgolubtsov/nonlinear-optimization-algorithms-multilang/master/nlp-unconstrained-cli/hooke-jeeves/js/src/woods.js|NLPUCCLIHooke.Woods}
 * @since    hooke-jeeves 0.1
 */
var Hooke = function() {

// === Public properties ======================================================

    /** Helper constants. */
    this.INDEX_ZERO      = 0;
    this.INDEX_ONE       = 1;
    var  ZERO_POINT_FIVE = 0.5; // <== Private.

    /**
     * Constant. The maximum number of variables.
     * @constant {Number} VARS
     * @public
     */
    this.VARS = 250;

    /**
     * Constant. The ending value of stepsize.
     * @constant {Number} EPSMIN
     * @public
     */
    this.EPSMIN = 1E-6;

    /**
     * Constant. The maximum number of iterations.
     * @constant {Number} IMAX
     * @public
     */
    this.IMAX = 5000;

    /**
     * The number of function evaluations.
     * @member  {Number} funEvals
     * @public
     * @default
     */
    this.funEvals = 0;

// === Private methods ========================================================

    /**
     * Helper method.
     * <br />
     * <br />Given a point, look for a better one nearby, one coord at a time.
     *
     * @function bestNearby
     *
     * @param {Number[]} delta    - The delta between prevBest and point.
     * @param {Number[]} point    - The coordinate from where to begin.
     * @param {Number}   prevBest - The previous best-valued coordinate.
     * @param {Number}   nVars    - The number of variables.
     * @param {Object}   f        - The user-supplied objective function
     *                              f(x,n).
     *
     * @returns {Number} The objective function value at a nearby.
     */
    var bestNearby = function(delta, point, prevBest, nVars, f) {
        var minF;
        var z = new Array(h.VARS);
        var fTmp;

        var i;

        minF = prevBest;

        for (i = 0; i < nVars; i++) {
            z[i] = point[i];
        }

        for (i = 0; i < nVars; i++) {
            z[i] = point[i] + delta[i];

            fTmp = f(z, nVars);

            if (fTmp < minF) {
                minF = fTmp;
            } else {
                delta[i] = 0.0 - delta[i];
                z[i]     = point[i] + delta[i];

                fTmp = f(z, nVars);

                if (fTmp < minF) {
                    minF = fTmp;
                } else {
                    z[i] = point[i];
                }
            }
        }

        for (i = 0; i < nVars; i++) {
            point[i] = z[i];
        }

        return minF;
    };

// === Public methods =========================================================

    /**
     * Main optimization method.
     * <br />
     * <br />The hooke subroutine itself.
     *
     * @function hooke
     * @public
     *
     * @param {Number}   nVars   - The number of variables.
     * @param {Number[]} startPt - The starting point coordinates.
     * @param {Number[]} endPt   - The ending point coordinates.
     * @param {Number}   rho     - The rho value.
     * @param {Number}   epsilon - The epsilon value.
     * @param {Number}   iterMax - The maximum number of iterations.
     * @param {Object}   f       - The user-supplied objective function f(x,n).
     *
     * @returns {Number} The number of iterations used
     *                   to find the local minimum.
     */
    this.hooke = function(nVars, startPt, endPt, rho, epsilon, iterMax, f) {
        var i;
        var iAdj;
        var iters;
        var j;
        var keep;

        var newX    = new Array(this.VARS);
        var xBefore = new Array(this.VARS);
        var delta   = new Array(this.VARS);
        var stepLength;
        var fBefore;
        var newF;
        var tmp;

        for (i = 0; i < nVars; i++) {
            xBefore[i] = startPt[i];
            newX[i]    = xBefore[i];

            delta[i] = Math.abs(startPt[i] * rho);

            if (delta[i] == 0.0) {
                delta[i] = rho;
            }
        }

        iAdj       = 0;
        stepLength = rho;
        iters      = 0;

        fBefore = f(newX, nVars);

        newF = fBefore;

        while ((iters < iterMax) && (stepLength > epsilon)) {
            iters++;
            iAdj++;

            console.log(
                "\nAfter "            + this.funEvals.toPrecision(4)
              + " funevals, f(x) =  " + fBefore.toExponential(4)
              + " at"
            );

            for (j = 0; j < nVars; j++) {
                console.log(
                    "   x[ " + j + "] = " + xBefore[j].toExponential(4)
                );
            }

            // Find best new point, one coord at a time.
            for (i = 0; i < nVars; i++) {
                newX[i] = xBefore[i];
            }

            newF = bestNearby(delta, newX, fBefore, nVars, f);

            // If we made some improvements, pursue that direction.
            keep = 1;

            while ((newF < fBefore) && (keep == 1)) {
                iAdj = 0;

                for (i = 0; i < nVars; i++) {
                    // Firstly, arrange the sign of delta[].
                    if (newX[i] <= xBefore[i]) {
                        delta[i] = 0.0 - Math.abs(delta[i]);
                    } else {
                        delta[i] = Math.abs(delta[i]);
                    }

                    // Now, move further in this direction.
                    tmp        = xBefore[i];
                    xBefore[i] = newX[i];
                    newX[i]    = newX[i] + newX[i] - tmp;
                }

                fBefore = newF;

                newF = bestNearby(delta, newX, fBefore, nVars, f);

                // If the further (optimistic) move was bad....
                if (newF >= fBefore) {
                    break;
                }

                /*
                 * Make sure that the differences between the new and the old
                 * points are due to actual displacements; beware of roundoff
                 * errors that might cause newF < fBefore.
                 */
                keep = 0;

                for (i = 0; i < nVars; i++) {
                    keep = 1;

                    if (Math.abs(newX[i] - xBefore[i])
                        > (ZERO_POINT_FIVE * Math.abs(delta[i]))) {

                        break;
                    } else {
                        keep = 0;
                    }
                }
            }

            if ((stepLength >= epsilon) && (newF >= fBefore)) {
                stepLength = stepLength * rho;

                for (i = 0; i < nVars; i++) {
                    delta[i] *= rho;
                }
            }
        }

        for (i = 0; i < nVars; i++) {
            endPt[i] = xBefore[i];
        }

        return iters;
    };
};

// Instantiating the Hooke class.
var h = NLPUCCLIHooke.Hooke = new Hooke();

// Exposing out public props and methods.
exports.INDEX_ZERO = h.INDEX_ZERO;
exports.INDEX_ONE  = h.INDEX_ONE;
exports.VARS       = h.VARS;
exports.EPSMIN     = h.EPSMIN;
exports.IMAX       = h.IMAX;
exports.funEvals   = h.funEvals;
exports.hooke      = h.hooke;

// vim:set nu et ts=4 sw=4:
