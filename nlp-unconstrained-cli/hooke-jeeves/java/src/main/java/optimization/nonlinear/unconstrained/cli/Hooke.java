/*
 * nlp-unconstrained-cli/hooke-jeeves/java/src/main/java/
 * optimization/nonlinear/unconstrained/cli/Hooke.java
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * Written by Radislav (Radicchio) Golubtsov, 2015-2023
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

package optimization.nonlinear.unconstrained.cli;

/**
 * The <code>Hooke</code> class contains methods for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 *
 * @author  Radislav (Radicchio) Golubtsov
 * @version 0.1.1
 * @see     optimization.nonlinear.unconstrained.cli.Rosenbrock
 * @see     optimization.nonlinear.unconstrained.cli.Woods
 * @since   hooke-jeeves 0.1
 */
public class Hooke {
    /** Constant. The maximum number of variables. */
    public static final int VARS = 250;

    /** Constant. The ending value of stepsize. */
    public static final double EPSMIN = 1E-6;

    /** Constant. The maximum number of iterations. */
    public static final int IMAX = 5000;

    /** Helper constants. */
    public  static final int    INDEX_ZERO      = 0;
    public  static final int    INDEX_ONE       = 1;
    private static final double ZERO_POINT_FIVE = 0.5;

    /** The number of function evaluations. */
    private static int funEvals = 0;

    /**
     * Getter for <code>funEvals</code>.
     *
     * @return The number of function evaluations.
     */
    public static int getFunEvals() {
        return funEvals;
    }

    /**
     * Setter for <code>funEvals</code>.
     *
     * @param __funEvals The number of function evaluations.
     */
    public static void setFunEvals(final int __funEvals) {
        funEvals = __funEvals;
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
    private double bestNearby(final double[] delta,
                              final double[] point,
                              final double prevBest,
                              final int nVars,
                              final Class objFunCls) {

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

            if (objFunCls.equals(Rosenbrock.class)) {
                fTmp = Rosenbrock.f(z, nVars);
            } else if (objFunCls.equals(Woods.class)) {
                fTmp = Woods.f(z, nVars);
            } else {
                fTmp = 0;
            }

            if (fTmp < minF) {
                minF = fTmp;
            } else {
                delta[i] = 0.0 - delta[i];
                z[i]     = point[i] + delta[i];

                if (objFunCls.equals(Rosenbrock.class)) {
                    fTmp = Rosenbrock.f(z, nVars);
                } else if (objFunCls.equals(Woods.class)) {
                    fTmp = Woods.f(z, nVars);
                } else {
                    fTmp = 0;
                }

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
     * Main optimization method.
     * <br />
     * <br />The hooke subroutine itself.
     *
     * @param nVars     The number of variables.
     * @param startPt   The starting point coordinates.
     * @param endPt     The ending point coordinates.
     * @param rho       The rho value.
     * @param epsilon   The epsilon value.
     * @param iterMax   The maximum number of iterations.
     * @param objFunCls The class in which the objective function is defined.
     *
     * @return The number of iterations used to find the local minimum.
     */
    public int hooke(final int nVars,
                     final double[] startPt,
                     final double[] endPt,
                     final double rho,
                     final double epsilon,
                     final int iterMax,
                     final Class objFunCls) {

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

        if (objFunCls.equals(Rosenbrock.class)) {
            fBefore = Rosenbrock.f(newX, nVars);
        } else if (objFunCls.equals(Woods.class)) {
            fBefore = Woods.f(newX, nVars);
        } else {
            fBefore = 0;
        }

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

            newF = bestNearby(delta, newX, fBefore, nVars, objFunCls);

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

                newF = bestNearby(delta, newX, fBefore, nVars, objFunCls);

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
    }

    /** Default constructor. */
    public Hooke() {}
}

// vim:set nu et ts=4 sw=4:
