/*
 * nlp-unconstrained-core/hooke-jeeves/java/hooke-jeeves/src/main/java/
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
 * Copyright (C) 2015 Radislav (Radic) Golubtsov
 */

/*
 * ============================================================================
   Nonlinear Optimization using the algorithm of Hooke and Jeeves
   12 February 1994 author: Mark G. Johnson

   Find a point X where the nonlinear function f(X) has a local  minimum.  X
   is an n-vector and f(X) is a scalar. In mathematical notation f:  R^n  ->
   R^1. The objective function f() is not required  to  be  continuous.  Nor
   does f() need to be differentiable. The program does not use  or  require
   derivatives of f().

   The software user supplies three things: a subroutine that computes f(X),
   an initial "starting guess" of the minimum point X, and  values  for  the
   algorithm convergence parameters. Then the program searches for  a  local
   minimum, beginning from the  starting  guess,  using  the  Direct  Search
   algorithm of Hooke and Jeeves.

   This C program is adapted from the Algol pseudocode found  in  "Algorithm
   178: Direct Search" by Arthur F. Kaupe Jr., Communications  of  the  ACM,
   Vol 6. p.313 (June 1963). It includes the improvements suggested by  Bell
   and Pike (CACM v.9, p.684, Sept 1966) and  those  of  Tomlin  and  Smith,
   "Remark on Algorithm 178" (CACM v.12). The original paper, which I  don't
   recommend as highly as the one by A.  Kaupe,  is:  R.  Hooke  and  T.  A.
   Jeeves, "Direct Search Solution of Numerical and  Statistical  Problems",
   Journal of the ACM, Vol.8, April 1961, pp.212-229.

   Calling sequence:
     int hooke(nvars, startpt, endpt, rho, epsilon, itermax)
       nvars    {an integer}           This is the number of  dimensions  in
                                       the domain of f(). It is  the  number
                                       of coordinates of the starting  point
                                       (and the minimum point.)
       startpt  {an array of doubles}  This is the  user-supplied  guess  at
                                       the minimum.
       endpt    {an array of doubles}  This is the  location  of  the  local
                                       minimum, calculated by the program
       rho      {a double}             This is a  user-supplied  convergence
                                       parameter (more detail below),  which
                                       should be set to a value between  0.0
                                       and 1.0. Larger values  of  rho  give
                                       greater probability of convergence on
                                       highly nonlinear functions, at a cost
                                       of more function evaluations. Smaller
                                       values of rho reduces the  number  of
                                       evaluations (and the program  running
                                       time),  but  increases  the  risk  of
                                       nonconvergence. See below.
       epsilon  {a double}             This is the criterion for halting the
                                       search  for  a  minimum.   When   the
                                       algorithm begins  to  make  less  and
                                       less progress on each  iteration,  it
                                       checks the halting criterion: if  the
                                       stepsize is below epsilon,  terminate
                                       the iteration and return the  current
                                       best estimate of the minimum.  Larger
                                       values of epsilon  (such  as  1.0e-4)
                                       give quicker running time, but a less
                                       accurate  estimate  of  the  minimum.
                                       Smaller values of  epsilon  (such  as
                                       1.0e-7) give longer running time, but
                                       a  more  accurate  estimate  of   the
                                       minimum.
       itermax  {an integer}           A  second,   rarely   used,   halting
                                       criterion. If the algorithm  uses  >=
                                       itermax iterations, halt.

   The user-supplied objective function f(x,n) should return a  C  "double".
   Its arguments are x -- an array of doubles, and n -- an integer. x is the
   point at which  f(x)  should  be  evaluated,  and  n  is  the  number  of
   coordinates of x. That is, n is the number of coefficients being fitted.

     rho, the algorithm convergence control

   The algorithm works by taking "steps" from one estimate of a minimum,  to
   another (hopefully better) estimate. Taking big steps gets to the minimum
   more quickly, at the risk of "stepping right over"  an  excellent  point.
   The stepsize is controlled by a user supplied parameter  called  rho.  At
   each iteration, the stepsize is multiplied by rho (0 < rho < 1),  so  the
   stepsize is successively reduced.
       Small values of rho correspond to big stepsize  changes,  which  make
   the algorithm run more quickly. However, there is  a  chance  (especially
   with highly nonlinear functions) that these big changes will accidentally
   overlook a promising search vector, leading to nonconvergence.
       Large values of rho correspond to small stepsize changes, which force
   the  algorithm  to   carefully   examine   nearby   points   instead   of
   optimistically  forging  ahead.  This   improves   the   probability   of
   convergence.
       The stepsize is reduced until  it  is  equal  to  (or  smaller  than)
   epsilon. So  the  number  of  iterations  performed  by  Hooke-Jeeves  is
   determined by rho and epsilon:

     rho ** (number_of_iterations) = epsilon

   In general it is a good idea to set rho to an  aggressively  small  value
   like 0.5 (hoping for fast convergence). Then, if the user  suspects  that
   the reported minimum is incorrect (or perhaps not accurate  enough),  the
   program can be run again with a larger value of rho such as  0.85,  using
   the result of the first minimization as the starting guess to  begin  the
   second minimization.

   Normal use:
     (1) Code your function f() in the C language
     (2) Install your starting guess {or read it in}
     (3) Run the program
     (4) {for the skeptical}: Use the computed minimum as the starting point
         for another run

   Data Fitting:
     Code your function f() to be the sum  of  the  squares  of  the  errors
     (differences) between the computed values and the measured values. Then
     minimize f() using Hooke-Jeeves.
         EXAMPLE: you have 20 datapoints (ti, yi) and you want to find A,B,C
     such that (A * t * t) + (B * exp(t)) + (C * tan(t)) fits  the  data  as
     closely as possible. Then f() is just

     f(x) = SUM(measured_y[i] - ((A * t[i] * t[i]) +
                                 (B * exp(t[i])) +
                                 (C * tan(t[i])))) ^ 2

     where x[] is a 3-vector consisting of {A, B, C}.

   The author of this software is M.G. Johnson.
   Permission to use, copy, modify, and distribute  this  software  for  any
   purpose without fee is hereby granted, provided that this  entire  notice
   is included in all copies of any software which is or includes a copy  or
   modification of this  software  and  in  all  copies  of  the  supporting
   documentation for such software. THIS SOFTWARE IS BEING PROVIDED "AS IS",
   WITHOUT ANY EXPRESS OR  IMPLIED  WARRANTY.  IN  PARTICULAR,  NEITHER  THE
   AUTHOR  NOR  AT&T  MAKE  ANY  REPRESENTATION  OR  WARRANTY  OF  ANY  KIND
   CONCERNING THE MERCHANTABILITY OF THIS SOFTWARE OR ITS  FITNESS  FOR  ANY
   PARTICULAR PURPOSE.
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
    private static final double MINUS_ONE_POINT_TWO = -1.2;

    /** Default constructor. */
    private Rosenbrock() { }

    /**
     * Rosenbrock's classic parabolic valley ("banana") function.
     *
     * @param x vector of variables
     * @param n number of variables
     *
     * @return objective function value
     */
    private static double f(final double[] x, final int n) {
        double a;
        double b;
        double c;

        funEvals++;
        a = x[INDEX_ZERO];
        b = x[INDEX_ONE];
        c = ONE_HUNDRED_POINT_ZERO * (b - (a * a)) * (b - (a * a));

        return c + ((ONE_POINT_ZERO - a) * (ONE_POINT_ZERO - a));
    }

    /**
     * Given a point, look for a better one nearby, one coord at a time.
     *
     * @param delta    delta coordinates
     * @param point    point coordinates
     * @param prevBest previous best value
     * @param nVars    number of variables
     *
     * @return best nearby value
     */
    private static double bestNearby(final double[] delta,
        final double[] point, final double prevBest, final int nVars) {

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
                z[i] = point[i] + delta[i];
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
     * @param nVars   number of variables
     * @param startPt starting point coordinates
     * @param endPt   ending point coordinates
     * @param rho     rho value
     * @param epsilon epsilon value
     * @param iterMax maximum number of iterations
     *
     * @return number of iterations actually spent
     */
    private static int hooke(final int nVars, final double[] startPt,
        final double[] endPt, final double rho, final double epsilon,
        final int iterMax) {

        int i;
        int iAdj;
        int iters;
        int j;
        int keep;
        double[] newX = new double[VARS];
        double[] xBefore = new double[VARS];
        double[] delta = new double[VARS];
        double stepLength;
        double fBefore;
        double newF;
        double tmp;

        for (i = 0; i < nVars; i++) {
            xBefore[i] = startPt[i];
            newX[i] = xBefore[i];
            delta[i] = Math.abs(startPt[i] * rho);

            if (delta[i] == 0.0) {
                delta[i] = rho;
            }
        }

        iAdj = 0;
        stepLength = rho;
        iters = 0;
        fBefore = f(newX, nVars);
        newF = fBefore;

        while ((iters < iterMax) && (stepLength > epsilon)) {
            iters++;
            iAdj++;
            System.out.printf("\nAfter %5d funevals, f(x) =  %.4e at\n",
                funEvals, fBefore);

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
                    tmp = xBefore[i];
                    xBefore[i] = newX[i];
                    newX[i] = newX[i] + newX[i] - tmp;
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
     * @param args array of command-line arguments
     */
    public static void main(final String[] args) {
        int nVars;
        int iterMax;
        int jj;
        int i;
        double[] startPt = new double[VARS];
        double rho;
        double epsilon;
        double[] endPt = new double[VARS];

        // Starting guess for Rosenbrock's test function.
        nVars = 2;
        startPt[INDEX_ZERO] = MINUS_ONE_POINT_TWO;
        startPt[INDEX_ONE] = ONE_POINT_ZERO;
        iterMax = IMAX;
        rho = RHO_BEGIN;
        epsilon = EPSMIN;
        jj = hooke(nVars, startPt, endPt, rho, epsilon, iterMax);
        System.out.printf("\n\n\nHOOKE USED %d ITERATIONS, AND RETURNED\n",
            jj);

        for (i = 0; i < nVars; i++) {
            System.out.printf("x[%3d] = %15.7e \n", i, endPt[i]);
        }
    }
}

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
