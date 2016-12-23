/*
 * nlp-unconstrained-cli/hooke-jeeves/js/src/woods.js
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

"use strict";

/**
 * The <code>NLPUCCLIHooke</code> namespace is used as a container
 * for the <code>Woods</code> class.
 *
 * @namespace NLPUCCLIHooke
 */
var NLPUCCLIHooke = NLPUCCLIHooke || {};

/**
 * The <code>Woods</code> class is responsible for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 * <br />
 * <br />The objective function in this case
 * is the so-called &quot;Woods&quot; function.
 *
 * @class    Woods
 * @memberof NLPUCCLIHooke
 * @author   Radislav (Radicchio) Golubtsov
 * @version  0.1
 * @see      {@link https://raw.githubusercontent.com/rgolubtsov/nonlinear-optimization-algorithms-multilang/master/nlp-unconstrained-cli/hooke-jeeves/js/src/hooke.js|NLPUCCLIHooke.Hooke}
 * @since    hooke-jeeves 0.1
 */
var Woods = function() {
    // Importing the "hooke.js" module (containing the Hooke class).
    var h = require("./hooke.js");

// === Private properties =====================================================

    /** Helper constants. */
    var INDEX_TWO   =  2;
    var INDEX_THREE =  3;
    var ONE_HUNDRED =  100;
    var NINETY      =  90;
    var TEN         =  10;
    var TEN_POINT   =  10.;
    var FOUR        =  4;
    var MINUS_THREE = -3;
    var MINUS_ONE   = -1;

    /**
     * Constant. The stepsize geometric shrink.
     * <br />
     * <br />The Hooke &amp; Jeeves algorithm works reasonably well
     * on Rosenbrock's function, but can fare worse on some standard
     * test functions, depending on rho. Here is an example that works well
     * when rho = 0.5, but fares poorly with rho = 0.6, and better again
     * with rho = 0.8.
     * @constant {Number} RHO_WOODS
     */
    var RHO_WOODS = 0.6;

// === Private methods ========================================================

    /**
     * The user-supplied objective function f(x,n).
     * <br />
     * <br />Woods &ndash; a la More, Garbow &amp; Hillstrom
     * (TOMS algorithm 566).
     *
     * @function f
     *
     * @param {Number[]} x - The point at which f(x) should be evaluated.
     * @param {Number}   n - The number of coordinates of x.
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

        h.funEvals++;

        s1 = x[h.INDEX_ONE] - x[h.INDEX_ZERO] * x[h.INDEX_ZERO];
        s2 = 1              - x[h.INDEX_ZERO];
        s3 = x[h.INDEX_ONE] - 1;

        t1 = x[INDEX_THREE] - x[INDEX_TWO]    * x[INDEX_TWO];
        t2 = 1              - x[INDEX_TWO];
        t3 = x[INDEX_THREE] - 1;

        t4 = s3 + t3;
        t5 = s3 - t3;

        return (ONE_HUNDRED * (s1 * s1) + s2 * s2
                   + NINETY * (t1 * t1) + t2 * t2
                      + TEN * (t4 * t4) + t5 * t5 / TEN_POINT);
    };

// === Public methods =========================================================

    /**
     * Main program function.
     *
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

        var startPt = new Array(h.VARS);
        var rho;
        var epsilon;
        var endPt   = new Array(h.VARS);

        // Starting guess test problem "Woods".
        nVars                 = FOUR;
        startPt[h.INDEX_ZERO] = MINUS_THREE;
        startPt[h.INDEX_ONE]  = MINUS_ONE;
        startPt[INDEX_TWO]    = MINUS_THREE;
        startPt[INDEX_THREE]  = MINUS_ONE;
        iterMax               = h.IMAX;
        rho                   = RHO_WOODS;
        epsilon               = h.EPSMIN;

        jj = h.hooke(nVars, startPt, endPt, rho, epsilon, iterMax, f);

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
NLPUCCLIHooke.Woods = new Woods();

// Firing up computations.
NLPUCCLIHooke.Woods.main();

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
