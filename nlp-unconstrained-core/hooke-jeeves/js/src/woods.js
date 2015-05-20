/*
 * nlp-unconstrained-core/hooke-jeeves/js/src/woods.js
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

/**
 * The <code>NLPUCCoreHooke</code> namespace is used as a "container"
 * for the <code>Woods</code> class.
 *
 * @namespace NLPUCCoreHooke
 */
var NLPUCCoreHooke = NLPUCCoreHooke || {};

/**
 * The <code>Woods</code> class is responsible for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 * <br />
 * <br />The objective function in this case
 * is so called "Woods" function.
 *
 * @class    Woods
 * @memberof NLPUCCoreHooke
 * @author   Radislav (Radic) Golubtsov
 * @version  0.1
 * @see      {@link https://github.com/rgolubtsov/nonlinear-optimization-algorithms-multilang/blob/master/nlp-unconstrained-core/hooke-jeeves/js/src/rosenbrock.js|NLPUCCoreHooke.Rosenbrock}
 * @since    0.1
 */
var Woods = function() {

// === Private properties =====================================================

    /**
     * The number of function evaluations.
     * @member  {Number} funEvals
     * @default
     */
    var funEvals = 0;

    /**
     * The maximum number of variables.
     * @constant {Number} VARS
     */
    var VARS = 250;

    /**
     * The stepsize geometric.
     * <br />
     * <br />The Hooke & Jeeves algorithm works reasonably well on Rosenbrock's
     * function, but can fare worse on some standard test functions,
     * depending on rho. Here is an example that works well when rho = 0.5,
     * but fares poorly with rho = 0.6, and better again with rho = 0.8.
     * @constant {Number} RHO_WOODS
     */
    var RHO_WOODS = 0.6;

    /**
     * The ending value of stepsize.
     * @constant {Number} EPSMIN
     */
    var EPSMIN = 1E-6;

    /**
     * The maximum number of iterations.
     * @constant {Number} IMAX
     */
    var IMAX = 5000;

    /**
     * Helper constant.
     * @constant {Number} INDEX_ZERO
     */
    var INDEX_ZERO = 0;

    /**
     * Helper constant.
     * @constant {Number} INDEX_ONE
     */
    var INDEX_ONE = 1;

    /**
     * Helper constant.
     * @constant {Number} INDEX_TWO
     */
    var INDEX_TWO = 2;

    /**
     * Helper constant.
     * @constant {Number} INDEX_THREE
     */
    var INDEX_THREE = 3;

    /**
     * Helper constant.
     * @constant {Number} ONE_HUNDRED
     */
    var ONE_HUNDRED = 100;

    /**
     * Helper constant.
     * @constant {Number} NINETY
     */
    var NINETY = 90;

    /**
     * Helper constant.
     * @constant {Number} TEN
     */
    var TEN = 10;

    /**
     * Helper constant.
     * @constant {Number} TEN_POINT
     */
    var TEN_POINT = 10.;

    /**
     * Helper constant.
     * @constant {Number} ZERO_POINT_FIVE
     */
    var ZERO_POINT_FIVE = 0.5;

    /**
     * Helper constant.
     * @constant {Number} FOUR
     */
    var FOUR = 4;

    /**
     * Helper constant.
     * @constant {Number} MINUS_THREE
     */
    var MINUS_THREE = -3;

    /**
     * Helper constant.
     * @constant {Number} MINUS_ONE
     */
    var MINUS_ONE = -1;

// === Private methods ========================================================

    /**
     * Woods -- a la More, Garbow & Hillstrom (TOMS algorithm 566).
     * @function f
     *
     * @param {Number[]} x - The vector of variables.
     * @param {Number}   n - The number of variables.
     *
     * @returns {Number} The objective function value.
     */
    var f = function(x, n) {
        var s1;
        var s2;
        var s3;
        var t1;
        var t2;
        var t3;
        var t4;
        var t5;

        funEvals++;

        s1 = x[INDEX_ONE] - x[INDEX_ZERO] * x[INDEX_ZERO];
        s2 = 1 - x[INDEX_ZERO];
        s3 = x[INDEX_ONE] - 1;
        t1 = x[INDEX_THREE] - x[INDEX_TWO] * x[INDEX_TWO];
        t2 = 1 - x[INDEX_TWO];
        t3 = x[INDEX_THREE] - 1;
        t4 = s3 + t3;
        t5 = s3 - t3;

        return (ONE_HUNDRED * (s1 * s1) + s2 * s2
                   + NINETY * (t1 * t1) + t2 * t2
                      + TEN * (t4 * t4) + t5 * t5 / TEN_POINT);
    };

    /**
     * Given a point, look for a better one nearby, one coord at a time.
     * @function bestNearby
     *
     * @param {Number[]} delta    - The delta coordinates.
     * @param {Number[]} point    - The point coordinates.
     * @param {Number}   prevBest - The previous best value.
     * @param {Number}   nVars    - The number of variables.
     *
     * @returns {Number} The best nearby value.
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
     * @function hooke
     *
     * @param {Number}   nVars   - The number of variables.
     * @param {Number[]} startPt - The starting point coordinates.
     * @param {Number[]} endPt   - The ending point coordinates.
     * @param {Number}   rho     - The rho value.
     * @param {Number}   epsilon - The epsilon value.
     * @param {Number}   iterMax - The maximum number of iterations.
     *
     * @returns {Number} The number of iterations actually spent.
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

// === Public methods =========================================================

    /**
     * Main program function.
     * @function main
     * @public
     *
     * @param {String[]} args - The array of command-line arguments.
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

        // Starting guess test problem "Woods".
        nVars                = FOUR;
        startPt[INDEX_ZERO]  = MINUS_THREE;
        startPt[INDEX_ONE]   = MINUS_ONE;
        startPt[INDEX_TWO]   = MINUS_THREE;
        startPt[INDEX_THREE] = MINUS_ONE;
        iterMax              = IMAX;
        rho                  = RHO_WOODS;
        epsilon              = EPSMIN;

        jj = hooke(nVars, startPt, endPt, rho, epsilon, iterMax);

        console.log("\n\n\nHOOKE USED " + jj + " ITERATIONS, AND RETURNED");

        for (i = 0; i < nVars; i++) {
            console.log(
                "x[  " + i + "] =   " + endPt[i].toExponential(7) + " "
            );
        }

        console.log("True answer: f(1, 1, 1, 1) = 0.");
    };
};

// Instantiating the Woods class.
NLPUCCoreHooke.Woods = new Woods();

// Firing up computations.
NLPUCCoreHooke.Woods.main();

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
