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
 * <br />
 * <br />The objective function in this case
 * is the Rosenbrock's parabolic valley function.
 *
 * @author  Radislav (Radic) Golubtsov
 * @version 0.1
 * @see     optimization.nonlinear.unconstrained.core.Hooke
 * @since   findMinimum-jeeves 0.1
 */
public final class Rosenbrock implements ObjectiveFunction{
    /** Constant. The stepsize geometric shrink. */
    private static final double RHO_BEGIN = 0.5;

    /** Helper constants. */
    private static final double ONE_HUNDRED_POINT_ZERO =  100.0;
    private static final double ONE_POINT_ZERO         =  1.0;
    private static final int    TWO                    =  2;
    private static final double MINUS_ONE_POINT_TWO    = -1.2;

    /**
     * The user-supplied objective function f(x,n).
     * <br />
     * <br />Represents here the Rosenbrock's classic parabolic valley
     * (&quot;banana&quot;) function.
     *
     * @param x The point at which f(x) should be evaluated.
     *
     * @return The objective function value.
     */
    public double findValueForArguments(final double[] x) {
        double a;
        double b;
        double c;

        a = x[Hooke.INDEX_ZERO];
        b = x[Hooke.INDEX_ONE];

        c = ONE_HUNDRED_POINT_ZERO * (b - (a * a)) * (b - (a * a));

        return (c + ((ONE_POINT_ZERO - a) * (ONE_POINT_ZERO - a)));
    }

    /**
     * Main program function.
     *
     * @param args The array of command-line arguments.
     */
    public static void main(final String[] args) {
        int nVars;
        int iterMax;
        int numberOfIterations;
        int i;

        double[] startPt = new double[Hooke.MAXIMUM_NUMBER_OF_VARIABLES];
        double rho;
        double epsilon;
        double[] endPt   = new double[Hooke.MAXIMUM_NUMBER_OF_VARIABLES];

        nVars                     = TWO;
        startPt[Hooke.INDEX_ZERO] = MINUS_ONE_POINT_TWO;
        startPt[Hooke.INDEX_ONE]  = ONE_POINT_ZERO;
        iterMax                   = Hooke.MAXIMUM_NUMBER_OF_ITERATIONS;
        rho                       = RHO_BEGIN;
        epsilon                   = Hooke.ENDING_VALUE_OF_STEPSIZE;

        Hooke hooke = new Hooke();

        numberOfIterations = hooke.findMinimum(
                nVars, startPt, endPt, rho, epsilon, iterMax, new Rosenbrock()
        );

        System.out.println(
            "\n\n\nHOOKE USED " + numberOfIterations + " ITERATIONS, AND RETURNED"
        );

        for (i = 0; i < nVars; i++) {
            System.out.printf("x[%3d] = %15.7e \n", i, endPt[i]);
        }
    }

}
