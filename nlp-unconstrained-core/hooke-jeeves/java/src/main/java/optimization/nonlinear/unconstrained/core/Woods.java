/*
 * nlp-unconstrained-core/hooke-jeeves/java/hooke-jeeves/src/main/java/
 * optimization/nonlinear/unconstrained/core/Woods.java
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * Copyright (C) 2015 Radislav (Radic) Golubtsov
 */

package optimization.nonlinear.unconstrained.core;

/**
 * The <code>Woods</code> class is responsible for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 * <p>
 * The objective function in this case is so called &quot;Woods&quot; function.
 *
 * @author  Radislav (Radic) Golubtsov
 * @version 0.1, 20150430
 * @see     optimization.nonlinear.unconstrained.core.Rosenbrock
 * @since   hooke-jeeves 0.1
 */
public final class Woods {
    /** Number of function evaluations. */
    private static int funEvals = 0;

    /** Constant: maximum number of variables. */
    private static final int VARS = 250;

    /**
     * Constant: stepsize geometric.
     * <p>
     * The Hooke &amp; Jeeves algorithm works reasonably well on Rosenbrock's
     * function, but can fare worse on some standard test functions,
     * depending on rho. Here is an example that works well when rho = 0.5,
     * but fares poorly with rho = 0.6, and better again with rho = 0.8.
     */
    private static final double RHO_WOODS = 0.6;

    /** Constant: ending value of stepsize. */
    private static final double EPSMIN = 1E-6;

    /** Constant: maximum number of iterations. */
    private static final int IMAX = 5000;

    /** Helper constant. */
    private static final int INDEX_ZERO = 0;

    /** Helper constant. */
    private static final int INDEX_ONE = 1;

    /** Helper constant. */
    private static final int INDEX_TWO = 2;

    /** Helper constant. */
    private static final int INDEX_THREE = 3;

    /** Helper constant. */
    private static final int ONE_HUNDRED = 100;

    /** Helper constant. */
    private static final int NINETY = 90;

    /** Helper constant. */
    private static final int TEN = 10;

    /** Helper constant. */
    private static final double TEN_POINT = 10.;

    /** Helper constant. */
    private static final double ZERO_POINT_FIVE = 0.5;

    /** Helper constant. */
    private static final int FOUR = 4;

    /** Helper constant. */
    private static final int MINUS_THREE = -3;

    /** Helper constant. */
    private static final int MINUS_ONE = -1;

    /** Default constructor. */
    private Woods() { }

    /**
     * Woods &ndash; a la More, Garbow &amp; Hillstrom (TOMS algorithm 566).
     *
     * @param x The vector of variables.
     * @param n The number of variables.
     *
     * @return The objective function value.
     */
    private static double f(final double[] x, final int n) {
        double s1;
        double s2;
        double s3;
        double t1;
        double t2;
        double t3;
        double t4;
        double t5;

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
    }

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
    private static double bestNearby(final double[] delta,
                                     final double[] point,
                                     final double prevBest,
                                     final int nVars) {

        double minF;
        double[] z = new double[VARS];
        double fTmp;

        int i;

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
    }

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
    private static int hooke(final int nVars,
                             final double[] startPt,
                             final double[] endPt,
                             final double rho,
                             final double epsilon,
                             final int iterMax) {

        int i;
        int iAdj;
        int iters;
        int j;
        int keep;

        double[] newX    = new double[VARS];
        double[] xBefore = new double[VARS];
        double[] delta   = new double[VARS];
        double stepLength;
        double fBefore;
        double newF;
        double tmp;

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

            System.out.printf(
                "\nAfter %5d funevals, f(x) =  %.4e at\n", funEvals, fBefore
            );

            for (j = 0; j < nVars; j++) {
                System.out.printf("   x[%2d] = %.4e\n", j, xBefore[j]);
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
    }

    /**
     * Main program function.
     *
     * @param args The array of command-line arguments.
     */
    public static void main(final String[] args) {
        int nVars;
        int iterMax;
        int jj;
        int i;

        double[] startPt = new double[VARS];
        double rho;
        double epsilon;
        double[] endPt   = new double[VARS];

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

        System.out.printf(
            "\n\n\nHOOKE USED %d ITERATIONS, AND RETURNED\n", jj
        );

        for (i = 0; i < nVars; i++) {
            System.out.printf("x[%3d] = %15.7e \n", i, endPt[i]);
        }

        System.out.println("True answer: f(1, 1, 1, 1) = 0.");
    }
}

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
