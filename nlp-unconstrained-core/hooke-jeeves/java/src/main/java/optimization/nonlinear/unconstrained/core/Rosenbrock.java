/*
 * nlp-unconstrained-core/hooke-jeeves/java/src/main/java/
 * optimization/nonlinear/unconstrained/core/Rosenbrock.java
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

package optimization.nonlinear.unconstrained.core;

/**
 * The <code>Rosenbrock</code> class is responsible for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 * <p>
 * The objective function in this case is the Rosenbrock's parabolic valley
 * function.
 *
 * @author  Radislav (Radic) Golubtsov
 * @version 0.1, 20150430
 * @see     optimization.nonlinear.unconstrained.core.Woods
 * @since   hooke-jeeves 0.1
 */
public final class Rosenbrock {
    /** Number of function evaluations. */
    private static int funEvals = 0;

    /** Constant: maximum number of variables. */
    private static final int VARS = 250;

    /** Constant: stepsize geometric. */
    private static final double RHO_BEGIN = 0.5;

    /** Constant: ending value of stepsize. */
    private static final double EPSMIN = 1E-6;

    /** Constant: maximum number of iterations. */
    private static final int IMAX = 5000;

    /** Helper constant. */
    private static final int INDEX_ZERO = 0;

    /** Helper constant. */
    private static final int INDEX_ONE = 1;

    /** Helper constant. */
    private static final double ONE_HUNDRED_POINT_ZERO = 100.0;

    /** Helper constant. */
    private static final double ONE_POINT_ZERO = 1.0;

    /** Helper constant. */
    private static final double ZERO_POINT_FIVE = 0.5;

    /** Helper constant. */
    private static final int TWO = 2;

    /** Helper constant. */
    private static final double MINUS_ONE_POINT_TWO = -1.2;

    /** Default constructor. */
    private Rosenbrock() { }

    /**
     * Rosenbrock's classic parabolic valley ("banana") function.
     *
     * @param x The vector of variables.
     * @param n The number of variables.
     *
     * @return The objective function value.
     */
    private static double f(final double[] x, final int n) {
        double a;
        double b;
        double c;

        funEvals++;

        a = x[INDEX_ZERO];
        b = x[INDEX_ONE];

        c = ONE_HUNDRED_POINT_ZERO * (b - (a * a)) * (b - (a * a));

        return (c + ((ONE_POINT_ZERO - a) * (ONE_POINT_ZERO - a)));
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

        // Starting guess for Rosenbrock's test function.
        nVars               = TWO;
        startPt[INDEX_ZERO] = MINUS_ONE_POINT_TWO;
        startPt[INDEX_ONE]  = ONE_POINT_ZERO;
        iterMax             = IMAX;
        rho                 = RHO_BEGIN;
        epsilon             = EPSMIN;

        jj = hooke(nVars, startPt, endPt, rho, epsilon, iterMax);

        System.out.printf(
            "\n\n\nHOOKE USED %d ITERATIONS, AND RETURNED\n", jj
        );

        for (i = 0; i < nVars; i++) {
            System.out.printf("x[%3d] = %15.7e \n", i, endPt[i]);
        }
    }
}

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
