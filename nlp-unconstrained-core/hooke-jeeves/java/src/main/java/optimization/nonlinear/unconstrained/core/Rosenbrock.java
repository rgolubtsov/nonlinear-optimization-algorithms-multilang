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
 * @since   hooke-jeeves 0.1
 */
public final class Rosenbrock {
    /** Constant. The stepsize geometric. */
    private static final double RHO_BEGIN = 0.5;

    /** Helper constant. */
    private static final double ONE_HUNDRED_POINT_ZERO = 100.0;

    /** Helper constant. */
    private static final double ONE_POINT_ZERO = 1.0;

    /** Helper constant. */
    private static final int TWO = 2;

    /** Helper constant. */
    private static final double MINUS_ONE_POINT_TWO = -1.2;

    /**
     * The user-supplied objective function f(x,n).
     * <br />
     * <br />Represents here the Rosenbrock's classic parabolic valley
     * (&quot;banana&quot;) function.
     *
     * @param x The point at which f(x) should be evaluated.
     * @param n The number of coordinates of x.
     *
     * @return The objective function value.
     */
    public static double f(final double[] x, final int n) {
        double a;
        double b;
        double c;

        Hooke.setFunEvals(Hooke.getFunEvals() + 1);

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
        int jj;
        int i;

        double[] startPt = new double[Hooke.VARS];
        double rho;
        double epsilon;
        double[] endPt   = new double[Hooke.VARS];

        // Starting guess for Rosenbrock's test function.
        nVars                     = TWO;
        startPt[Hooke.INDEX_ZERO] = MINUS_ONE_POINT_TWO;
        startPt[Hooke.INDEX_ONE]  = ONE_POINT_ZERO;
        iterMax                   = Hooke.IMAX;
        rho                       = RHO_BEGIN;
        epsilon                   = Hooke.EPSMIN;

        // Instantiating the Hooke class.
        Hooke h = new Hooke();

        jj = h.hooke(
            nVars, startPt, endPt, rho, epsilon, iterMax, new Rosenbrock()
        );

        System.out.println(
            "\n\n\nHOOKE USED " + jj + " ITERATIONS, AND RETURNED"
        );

        for (i = 0; i < nVars; i++) {
            System.out.printf("x[%3d] = %15.7e \n", i, endPt[i]);
        }
    }

    /** Default constructor. */
    public Rosenbrock() {}
}

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
