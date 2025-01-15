/*
 * nlp-unconstrained-cli/hooke-jeeves/vala/src/hooke.vala
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * Written by Radislav (Radicchio) Golubtsov, 2015-2025
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

namespace CLIHooke {

/** Constant: The maximum number of variables. */
const uint VARS = 250;

/** Constant: The stepsize geometric shrink. */
const double RHO_BEGIN = 0.5;

/**
 * Constant: The stepsize geometric shrink.
 * The Hooke &amp; Jeeves algorithm works reasonably well
 * on Rosenbrock's function, but can fare worse on some standard
 * test functions, depending on rho. Here is an example that works well
 * when rho = 0.5, but fares poorly with rho = 0.6, and better again
 * with rho = 0.8.
 */
const double RHO_WOODS = 0.6;

/** Constant: The ending value of stepsize. */
const double EPSMIN = 1E-6;

/** Constant: The maximum number of iterations. */
const uint IMAX = 5000;

// Helper constants.
const uint   INDEX_ZERO          =  0;
const uint   INDEX_ONE           =  1;
const uint   INDEX_TWO           =  2;
const uint   INDEX_THREE         =  3;
const uint   TWO                 =  2;
const uint   FOUR                =  4;
const double MINUS_ONE_POINT_TWO = -1.2;
const double ONE_POINT_ZERO      =  1.0;
const int    MINUS_THREE         = -3;
const int    MINUS_ONE           = -1;
const double ZERO_POINT_FIVE     =  0.5;

/**
 * Helper function.
 * Given a point, look for a better one nearby, one coord at a time.
 *
 * @param delta    The delta between <code>prevbest</code>
 *                 and <code>point</code>.
 * @param point    The coordinate from where to begin.
 * @param prevbest The previous best-valued coordinate.
 * @param nvars    The number of variables.
 * @param funevals The number of function evaluations container (FunEvals *).
 *
 * @return The objective function value at a nearby.
 */
double best_nearby(double   *delta,
                   double   *point,
                   double    prevbest,
                   uint      nvars,
                   FunEvals *funevals) {

    double minf;

//==double z[VARS]; // <== Just to remember: Vala 0.22 failed on this,
    double z[250];  // .oO <== but succeeded with this one.))

    double ftmp;

    uint i;

    minf = prevbest;

    for (i = 0; i < nvars; i++) {
        z[i] = point[i];
    }

    for (i = 0; i < nvars; i++) {
        z[i] = point[i] + delta[i];

        ftmp = f(z, nvars, funevals);

        if (ftmp < minf) {
            minf = ftmp;
        } else {
            delta[i] = 0.0 - delta[i];
            z[i]     = point[i] + delta[i];

            ftmp = f(z, nvars, funevals);

            if (ftmp < minf) {
                minf = ftmp;
            } else {
                z[i] = point[i];
            }
        }
    }

    for (i = 0; i < nvars; i++) {
        point[i] = z[i];
    }

    return minf;
}

/**
 * Main optimization function.
 * The hooke subroutine itself.
 *
 * @param nvars   The number of variables.
 * @param startpt The starting point coordinates.
 * @param endpt   The ending point coordinates.
 * @param rho     The rho value.
 * @param epsilon The epsilon value.
 * @param itermax The maximum number of iterations.
 *
 * @return The number of iterations used to find the local minimum.
 */
uint hooke(uint    nvars,
           double *startpt,
           double *endpt,
           double  rho,
           double  epsilon,
           uint    itermax) {

    uint i;
    uint iadj;
    uint iters;
    uint j;
    uint keep;

//==double newx[VARS]; // <== Just to remember: Vala 0.22 failed on this,
    double newx[250];  // .oO <== but succeeded with this one.))

//==double xbefore[VARS]; // <== Just to remember: Vala 0.22 failed on this,
    double xbefore[250];  // .oO <== but succeeded with this one.))

//==double delta[VARS]; // <== Just to remember: Vala 0.22 failed on this,
    double delta[250];  // .oO <== but succeeded with this one.))

    double steplength;
    double fbefore;
    double newf;
    double tmp;

    for (i = 0; i < nvars; i++) {
        newx[i] = xbefore[i] = startpt[i];

        delta[i] = Math.fabs(startpt[i] * rho);

        if (delta[i] == 0.0) {
            delta[i] = rho;
        }
    }

    iadj       = 0;
    steplength = rho;
    iters      = 0;

    // Instantiating the FunEvals class.
    FunEvals *funevals = new FunEvals();

    fbefore = f(newx, nvars, funevals);

    newf = fbefore;

    while ((iters < itermax) && (steplength > epsilon)) {
        iters++;
        iadj++;

        stdout.printf("\nAfter %5u funevals, f(x) =  %.4le at\n",
            funevals->get_funevals(), fbefore);

        for (j = 0; j < nvars; j++) {
            stdout.printf("   x[%2u] = %.4le\n", j, xbefore[j]);
        }

        // Find best new point, one coord at a time.
        for (i = 0; i < nvars; i++) {
            newx[i] = xbefore[i];
        }

        newf = best_nearby(delta, newx, fbefore, nvars, funevals);

        // If we made some improvements, pursue that direction.
        keep = 1;

        while ((newf < fbefore) && (keep == 1)) {
            iadj = 0;

            for (i = 0; i < nvars; i++) {
                // Firstly, arrange the sign of delta[].
                if (newx[i] <= xbefore[i]) {
                    delta[i] = 0.0 - Math.fabs(delta[i]);
                } else {
                    delta[i] = Math.fabs(delta[i]);
                }

                // Now, move further in this direction.
                tmp        = xbefore[i];
                xbefore[i] = newx[i];
                newx[i]    = newx[i] + newx[i] - tmp;
            }

            fbefore = newf;

            newf = best_nearby(delta, newx, fbefore, nvars, funevals);

            // If the further (optimistic) move was bad....
            if (newf >= fbefore) {
                break;
            }

            /*
             * Make sure that the differences between the new and the old
             * points are due to actual displacements; beware of roundoff
             * errors that might cause newf < fbefore.
             */
            keep = 0;

            for (i = 0; i < nvars; i++) {
                keep = 1;

                if (Math.fabs(newx[i] - xbefore[i])
                    > (ZERO_POINT_FIVE * Math.fabs(delta[i]))) {

                    break;
                } else {
                    keep = 0;
                }
            }
        }

        if ((steplength >= epsilon) && (newf >= fbefore)) {
            steplength = steplength * rho;

            for (i = 0; i < nvars; i++) {
                delta[i] *= rho;
            }
        }
    }

    for (i = 0; i < nvars; i++) {
        endpt[i] = xbefore[i];
    }

    // Destroying the FunEvals class instance.
    delete funevals;

    return iters;
}

} // namespace CLIHooke

// Main program function main() :-).
public static int main() {
    uint nvars;
    uint itermax;
    uint jj;
    uint i;

//==double startpt[CLIHooke.VARS]; // <==Just to rmmbr:Vala0.22 failed on this,
    double startpt[250];           // .oO <== but succeeded with this one.))

    double rho;
    double epsilon;

//==double endpt[CLIHooke.VARS]; // <==Just to rmmbr: Vala 0.22 failed on this,
    double endpt[250];           // .oO <== but succeeded with this one.))

#if (!WOODS)
    // Starting guess for Rosenbrock's test function.
    nvars                         = CLIHooke.TWO;
    startpt[CLIHooke.INDEX_ZERO]  = CLIHooke.MINUS_ONE_POINT_TWO;
    startpt[CLIHooke.INDEX_ONE]   = CLIHooke.ONE_POINT_ZERO;
    rho                           = CLIHooke.RHO_BEGIN;
#else
    // Starting guess test problem "Woods".
    nvars                         = CLIHooke.FOUR;
    startpt[CLIHooke.INDEX_ZERO]  = CLIHooke.MINUS_THREE;
    startpt[CLIHooke.INDEX_ONE]   = CLIHooke.MINUS_ONE;
    startpt[CLIHooke.INDEX_TWO]   = CLIHooke.MINUS_THREE;
    startpt[CLIHooke.INDEX_THREE] = CLIHooke.MINUS_ONE;
    rho                           = CLIHooke.RHO_WOODS;
#endif

    itermax = CLIHooke.IMAX;
    epsilon = CLIHooke.EPSMIN;

    jj = CLIHooke.hooke(nvars, startpt, endpt, rho, epsilon, itermax);

    stdout.printf("\n\n\nHOOKE USED %u ITERATIONS, AND RETURNED\n", jj);

    for (i = 0; i < nvars; i++) {
        stdout.printf("x[%3u] = %15.7le \n", i, endpt[i]);
    }

#if (WOODS)
    stdout.puts("True answer: f(1, 1, 1, 1) = 0.\n");
#endif

    return Posix.EXIT_SUCCESS;
}

// vim:set nu et ts=4 sw=4:
