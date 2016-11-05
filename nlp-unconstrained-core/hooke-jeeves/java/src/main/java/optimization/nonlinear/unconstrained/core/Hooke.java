/*
 * nlp-unconstrained-core/hooke-jeeves/java/src/main/java/
 * optimization/nonlinear/unconstrained/core/Hooke.java
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
 * The <code>Hooke</code> class contains methods for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 *
 * @author  Radislav (Radic) Golubtsov
 * @version 0.1
 * @see     optimization.nonlinear.unconstrained.core.Rosenbrock
 * @see     optimization.nonlinear.unconstrained.core.Woods
 * @since   findMinimum-jeeves 0.1
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
                              final ObjectiveFunction objFunCls) {

        double minF;
        double[] z = new double[VARS];
        double currentFunctionValue;

        int iterationNumber;

        minF = prevBest;

        for (iterationNumber = 0; iterationNumber < nVars; iterationNumber++) {
            z[iterationNumber] = point[iterationNumber];
        }

        for (iterationNumber = 0; iterationNumber < nVars; iterationNumber++) {
            z[iterationNumber] = point[iterationNumber] + delta[iterationNumber];
            currentFunctionValue = objFunCls.objectiveFunctionValue(z);
//            is this necessary ? Maybe use Java 8 default method returning 0 then?
//        } else {
//            fBefore = 0;
//        }

            if (currentFunctionValue < minF) {
                minF = currentFunctionValue;
            } else {
                delta[iterationNumber] = 0.0 - delta[iterationNumber];
                z[iterationNumber]     = point[iterationNumber] + delta[iterationNumber];
                currentFunctionValue = objFunCls.objectiveFunctionValue(z);
//            is this necessary ? Maybe use Java 8 default method returning 0 then?
//        } else {
//            fBefore = 0;
//        }

                if (currentFunctionValue < minF) {
                    minF = currentFunctionValue;
                } else {
                    z[iterationNumber] = point[iterationNumber];
                }
            }
        }

        for (iterationNumber = 0; iterationNumber < nVars; iterationNumber++) {
            point[iterationNumber] = z[iterationNumber];
        }

        return minF;
    }

    /**
     * Main optimization method.
     * <br />
     * <br />The findMinimum subroutine itself.
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
    public int findMinimum(int nVars,
                           double[] startPt,
                           double[] endPt,
                           double rho,
                           double epsilon,
                           int iterMax,
                           ObjectiveFunction objFunCls) {

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

        fBefore=objFunCls.objectiveFunctionValue(newX);
//            is this necessary ? Maybe use Java 8 default method returning 0 then?
//        } else {
//            fBefore = 0;
//        }

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

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
