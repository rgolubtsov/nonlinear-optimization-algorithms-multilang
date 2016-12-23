/*
 * nlp-unconstrained-cli/hooke-jeeves/js/src/rosenbrock.js
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
 * for the <code>Rosenbrock</code> class.
 *
 * @namespace NLPUCCLIHooke
 */
var NLPUCCLIHooke = NLPUCCLIHooke || {};

/**
 * The <code>Rosenbrock</code> class is responsible for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 * <br />
 * <br />The objective function in this case
 * is the Rosenbrock's parabolic valley function.
 *
 * @class    Rosenbrock
 * @memberof NLPUCCLIHooke
 * @author   Radislav (Radicchio) Golubtsov
 * @version  0.1
 * @see      {@link https://raw.githubusercontent.com/rgolubtsov/nonlinear-optimization-algorithms-multilang/master/nlp-unconstrained-cli/hooke-jeeves/js/src/hooke.js|NLPUCCLIHooke.Hooke}
 * @since    hooke-jeeves 0.1
 */
var Rosenbrock = function() {
    // Importing the "hooke.js" module (containing the Hooke class).
    var h = require("./hooke.js");

// === Private properties =====================================================

    /** Helper constants. */
    var ONE_HUNDRED_POINT_ZERO =  100.0;
    var ONE_POINT_ZERO         =  1.0;
    var TWO                    =  2;
    var MINUS_ONE_POINT_TWO    = -1.2;

    /**
     * Constant. The stepsize geometric shrink.
     * @constant {Number} RHO_BEGIN
     */
    var RHO_BEGIN = 0.5;

// === Private methods ========================================================

    /**
     * The user-supplied objective function f(x,n).
     * <br />
     * <br />Represents here the Rosenbrock's classic parabolic valley
     * (&quot;banana&quot;) function.
     *
     * @function f
     *
     * @param {Number[]} x - The point at which f(x) should be evaluated.
     * @param {Number}   n - The number of coordinates of x.
     *
     * @returns {Number} The objective function value.
     */
    var f = function(x, n) {
        var a;
        var b;
        var c;

        h.funEvals++;

        a = x[h.INDEX_ZERO];
        b = x[h.INDEX_ONE];

        c = ONE_HUNDRED_POINT_ZERO * (b - (a * a)) * (b - (a * a));

        return (c + ((ONE_POINT_ZERO - a) * (ONE_POINT_ZERO - a)));
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

        // Starting guess for Rosenbrock's test function.
        nVars                 = TWO;
        startPt[h.INDEX_ZERO] = MINUS_ONE_POINT_TWO;
        startPt[h.INDEX_ONE]  = ONE_POINT_ZERO;
        iterMax               = h.IMAX;
        rho                   = RHO_BEGIN;
        epsilon               = h.EPSMIN;

        jj = h.hooke(nVars, startPt, endPt, rho, epsilon, iterMax, f);

        console.log("\n\n\nHOOKE USED " + jj + " ITERATIONS, AND RETURNED");

        for (i = 0; i < nVars; i++) {
            console.log(
                "x[  " + i + "] =   " + endPt[i].toExponential(7) + " "
            );
        }
    };
};

// Instantiating the Rosenbrock class.
NLPUCCLIHooke.Rosenbrock = new Rosenbrock();

// Firing up computations.
NLPUCCLIHooke.Rosenbrock.main();

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
