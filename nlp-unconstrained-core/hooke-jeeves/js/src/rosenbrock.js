/*
 * nlp-unconstrained-core/hooke-jeeves/js/src/rosenbrock.js
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * Written by Radislav (Radic) Golubtsov <https://github.com/rgolubtsov>
 */

/** Default namespace. */
var NLPUCCoreHooke = NLPUCCoreHooke || {};

/**
 * The Rosenbrock class is responsible for solving a nonlinear optimization
 * problem using the algorithm of Hooke and Jeeves.
 *
 * The objective function in this case is the Rosenbrock's parabolic valley
 * function.
 *
 * @constructor
 * @author  Radislav (Radic) Golubtsov
 * @version 0.1
 * @see     NLPUCCoreHooke.Woods
 */
var Rosenbrock = function() {
    /** Number of function evaluations. */
    var funEvals = 0;

    /** Constant: maximum number of variables. */
    var VARS = 250;

    /** Constant: stepsize geometric. */
    var RHO_BEGIN = 0.5;

    /** Constant: ending value of stepsize. */
    var EPSMIN = 1E-6;

    /** Constant: maximum number of iterations. */
    var IMAX = 5000;

    /** Helper constant. */
    var INDEX_ZERO = 0;

    /** Helper constant. */
    var INDEX_ONE = 1;

    /** Helper constant. */
    var ONE_HUNDRED_POINT_ZERO = 100.0;

    /** Helper constant. */
    var ONE_POINT_ZERO = 1.0;

    /** Helper constant. */
    var ZERO_POINT_FIVE = 0.5;

    /** Helper constant. */
    var TWO = 2;

    /** Helper constant. */
    var MINUS_ONE_POINT_TWO = -1.2;

    /**
     * Rosenbrock's classic parabolic valley ("banana") function.
     *
     * @param x The vector of variables.
     * @param n The number of variables.
     *
     * @return The objective function value.
     */
    var f = function(x, n) {
        var a;
        var b;
        var c;

        funEvals++;

        a = x[INDEX_ZERO];
        b = x[INDEX_ONE];

        c = ONE_HUNDRED_POINT_ZERO * (b - (a * a)) * (b - (a * a));

        return (c + ((ONE_POINT_ZERO - a) * (ONE_POINT_ZERO - a)));
    };

    /**
     * Given a point, look for a better one nearby, one coord at a time.
     *
     * @param delta    The delta coordinates.
     * @param point    The point coordinates.
     * @param prevBest The previous best value.
     * @param nVars    The number of variables.
     *
     * @return The best nearby value.
     */
    var bestNearby = function(delta, point, prevBest, nVars) {
        var minF;
        var z = new Array(VARS);
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

    /**
     * The hooke subroutine itself.
     *
     * @param nVars   The number of variables.
     * @param startPt The starting point coordinates.
     * @param endPt   The ending point coordinates.
     * @param rho     The rho value.
     * @param epsilon The epsilon value.
     * @param iterMax The maximum number of iterations.
     *
     * @return The number of iterations actually spent.
     */
    var hooke = function(nVars, startPt, endPt, rho, epsilon, iterMax) {
        var i;
        var iAdj;
        var iters;
        var j;
        var keep;

        var newX    = new Array(VARS);
        var xBefore = new Array(VARS);
        var delta   = new Array(VARS);
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
                "\nAfter "            + funEvals.toPrecision(4)
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

            newF = bestNearby(delta, newX, fBefore, nVars);

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

                newF = bestNearby(delta, newX, fBefore, nVars);

                // If the further (optimistic) move was bad....
                if (newF >= fBefore) {
                    break;
                }

                /*
                 * Make sure that the differences between the new and the old
                 * points are due to actual displacements; beware of roundoff
                 * errors that might cause newf < fbefore.
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

    /**
     * Main program function.
     *
     * @param args The array of command-line arguments.
     */
    this.main = function(args) {
        var nVars;
        var iterMax;
        var jj;
        var i;

        var startPt = new Array(VARS);
        var rho;
        var epsilon;
        var endPt   = new Array(VARS);

        // Starting guess for Rosenbrock's test function.
        nVars               = TWO;
        startPt[INDEX_ZERO] = MINUS_ONE_POINT_TWO;
        startPt[INDEX_ONE]  = ONE_POINT_ZERO;
        iterMax             = IMAX;
        rho                 = RHO_BEGIN;
        epsilon             = EPSMIN;

        jj = hooke(nVars, startPt, endPt, rho, epsilon, iterMax);

        console.log("\n\n\nHOOKE USED " + jj + " ITERATIONS, AND RETURNED");

        for (i = 0; i < nVars; i++) {
            console.log(
                "x[  " + i + "] =   " + endPt[i].toExponential(7) + " "
            );
        }
    };
};

// Instantiating the Rosenbrock class.
NLPUCCoreHooke.Rosenbrock = new Rosenbrock();

// Firing up computations.
NLPUCCoreHooke.Rosenbrock.main();

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
