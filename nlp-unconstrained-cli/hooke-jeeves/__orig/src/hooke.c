/*
 * nlp-unconstrained-cli/hooke-jeeves/__orig/src/hooke.c
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * (Re-)Written by Radislav (Radicchio) Golubtsov, 2016
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

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Nonlinear Optimization using the algorithm of Hooke and Jeeves            *
 * 12 February 1994 author: Mark G. Johnson                                  *
 *                                                                           *
 * Find a point X where the nonlinear function f(X) has a local  minimum.  X *
 * is an n-vector and f(X) is a scalar. In mathematical notation f:  R^n  -> *
 * R^1. The objective function f() is not required  to  be  continuous.  Nor *
 * does f() need to be differentiable. The program does not use  or  require *
 * derivatives of f().                                                       *
 *                                                                           *
 * The software user supplies three things: a subroutine that computes f(X), *
 * an initial "starting guess" of the minimum point X, and  values  for  the *
 * algorithm convergence parameters. Then the program searches for  a  local *
 * minimum, beginning from the  starting  guess,  using  the  Direct  Search *
 * algorithm of Hooke and Jeeves.                                            *
 *                                                                           *
 * This C program is adapted from the Algol pseudocode found  in  "Algorithm *
 * 178: Direct Search" by Arthur F. Kaupe Jr., Communications  of  the  ACM, *
 * Vol 6. p.313 (June 1963). It includes the improvements suggested by  Bell *
 * and Pike (CACM v.9, p.684, Sept 1966) and  those  of  Tomlin  and  Smith, *
 * "Remark on Algorithm 178" (CACM v.12). The original paper, which I  don't *
 * recommend as highly as the one by A.  Kaupe,  is:  R.  Hooke  and  T.  A. *
 * Jeeves, "Direct Search Solution of Numerical and  Statistical  Problems", *
 * Journal of the ACM, Vol.8, April 1961, pp.212-229.                        *
 *                                                                           *
 * Calling sequence:                                                         *
 *   int hooke(nvars, startpt, endpt, rho, epsilon, itermax)                 *
 *     nvars    {an integer}           This is the number of  dimensions  in *
 *                                     the domain of f(). It is  the  number *
 *                                     of coordinates of the starting  point *
 *                                     (and the minimum point.)              *
 *     startpt  {an array of doubles}  This is the  user-supplied  guess  at *
 *                                     the minimum.                          *
 *     endpt    {an array of doubles}  This is the  location  of  the  local *
 *                                     minimum, calculated by the program    *
 *     rho      {a double}             This is a  user-supplied  convergence *
 *                                     parameter (more detail below),  which *
 *                                     should be set to a value between  0.0 *
 *                                     and 1.0. Larger values  of  rho  give *
 *                                     greater probability of convergence on *
 *                                     highly nonlinear functions, at a cost *
 *                                     of more function evaluations. Smaller *
 *                                     values of rho reduces the  number  of *
 *                                     evaluations (and the program  running *
 *                                     time),  but  increases  the  risk  of *
 *                                     nonconvergence. See below.            *
 *     epsilon  {a double}             This is the criterion for halting the *
 *                                     search  for  a  minimum.   When   the *
 *                                     algorithm begins  to  make  less  and *
 *                                     less progress on each  iteration,  it *
 *                                     checks the halting criterion: if  the *
 *                                     stepsize is below epsilon,  terminate *
 *                                     the iteration and return the  current *
 *                                     best estimate of the minimum.  Larger *
 *                                     values of epsilon  (such  as  1.0e-4) *
 *                                     give quicker running time, but a less *
 *                                     accurate  estimate  of  the  minimum. *
 *                                     Smaller values of  epsilon  (such  as *
 *                                     1.0e-7) give longer running time, but *
 *                                     a  more  accurate  estimate  of   the *
 *                                     minimum.                              *
 *     itermax  {an integer}           A  second,   rarely   used,   halting *
 *                                     criterion. If the algorithm  uses  >= *
 *                                     itermax iterations, halt.             *
 *                                                                           *
 * The user-supplied objective function f(x,n) should return a  C  "double". *
 * Its arguments are x -- an array of doubles, and n -- an integer. x is the *
 * point at which  f(x)  should  be  evaluated,  and  n  is  the  number  of *
 * coordinates of x. That is, n is the number of coefficients being fitted.  *
 *                                                                           *
 *   rho, the algorithm convergence control                                  *
 *                                                                           *
 * The algorithm works by taking "steps" from one estimate of a minimum,  to *
 * another (hopefully better) estimate. Taking big steps gets to the minimum *
 * more quickly, at the risk of "stepping right over"  an  excellent  point. *
 * The stepsize is controlled by a user supplied parameter  called  rho.  At *
 * each iteration, the stepsize is multiplied by rho (0 < rho < 1),  so  the *
 * stepsize is successively reduced.                                         *
 *     Small values of rho correspond to big stepsize  changes,  which  make *
 * the algorithm run more quickly. However, there is  a  chance  (especially *
 * with highly nonlinear functions) that these big changes will accidentally *
 * overlook a promising search vector, leading to nonconvergence.            *
 *     Large values of rho correspond to small stepsize changes, which force *
 * the  algorithm  to   carefully   examine   nearby   points   instead   of *
 * optimistically  forging  ahead.  This   improves   the   probability   of *
 * convergence.                                                              *
 *     The stepsize is reduced until  it  is  equal  to  (or  smaller  than) *
 * epsilon. So  the  number  of  iterations  performed  by  Hooke-Jeeves  is *
 * determined by rho and epsilon:                                            *
 *                                                                           *
 *   rho ** (number_of_iterations) = epsilon                                 *
 *                                                                           *
 * In general it is a good idea to set rho to an  aggressively  small  value *
 * like 0.5 (hoping for fast convergence). Then, if the user  suspects  that *
 * the reported minimum is incorrect (or perhaps not accurate  enough),  the *
 * program can be run again with a larger value of rho such as  0.85,  using *
 * the result of the first minimization as the starting guess to  begin  the *
 * second minimization.                                                      *
 *                                                                           *
 * Normal use:                                                               *
 *   (1) Code your function f() in the C language                            *
 *   (2) Install your starting guess {or read it in}                         *
 *   (3) Run the program                                                     *
 *   (4) {for the skeptical}: Use the computed minimum as the starting point *
 *       for another run                                                     *
 *                                                                           *
 * Data Fitting:                                                             *
 *   Code your function f() to be the sum  of  the  squares  of  the  errors *
 *   (differences) between the computed values and the measured values. Then *
 *   minimize f() using Hooke-Jeeves.                                        *
 *   ==> EXAMPLE: you have 20 datapoints (ti, yi) and you want to find A,B,C *
 *   such that (A * t * t) + (B * exp(t)) + (C * tan(t)) fits  the  data  as *
 *   closely as possible. Then f() is just                                   *
 *                                                                           *
 *   f(x) = SUM(measured_y[i] - ((A * t[i] * t[i]) +                         *
 *                               (B * exp(t[i]))   +                         *
 *                               (C * tan(t[i])))) ^ 2                       *
 *                                                                           *
 *   where x[] is a 3-vector consisting of {A, B, C}.                        *
 *                                                                           *
 * The author of this software is M.G. Johnson.                              *
 * Permission to use, copy, modify, and distribute  this  software  for  any *
 * purpose without fee is hereby granted, provided that this  entire  notice *
 * is included in all copies of any software which is or includes a copy  or *
 * modification of this  software  and  in  all  copies  of  the  supporting *
 * documentation for such software. THIS SOFTWARE IS BEING PROVIDED "AS IS", *
 * WITHOUT ANY EXPRESS OR  IMPLIED  WARRANTY.  IN  PARTICULAR,  NEITHER  THE *
 * AUTHOR  NOR  AT&T  MAKE  ANY  REPRESENTATION  OR  WARRANTY  OF  ANY  KIND *
 * CONCERNING THE MERCHANTABILITY OF THIS SOFTWARE OR ITS  FITNESS  FOR  ANY *
 * PARTICULAR PURPOSE.                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <math.h>

/** Constant. The maximum number of variables. */
#define VARS 250

#ifndef WOODS
/** Constant. The stepsize geometric shrink. */
#define RHO_BEGIN 0.5
#else
/**
 * Constant. The stepsize geometric shrink.
 * <br />
 * <br />The Hooke &amp; Jeeves algorithm works reasonably well
 * on Rosenbrock's function, but can fare worse on some standard
 * test functions, depending on rho. Here is an example that works well
 * when rho = 0.5, but fares poorly with rho = 0.6, and better again
 * with rho = 0.8.
 */
#define RHO_WOODS 0.6
#endif

/** Constant. The ending value of stepsize. */
#define EPSMIN 1E-6

/** Constant. The maximum number of iterations. */
#define IMAX 5000

/** The number of function evaluations. */
int FUNEVALS = 0;

/* Forward declaration for the objective function. */
double
f();

/**
 * Helper function.
 * <br />
 * <br />Given a point, look for a better one nearby, one coord at a time.
 *
 * @param delta    The delta between <code>prevbest</code>
 *                 and <code>point</code>.
 * @param point    The coordinate from where to begin.
 * @param prevbest The previous best-valued coordinate.
 * @param nvars    The number of variables.
 *
 * @return The objective function value at a nearby.
 */
double
best_nearby(delta, point, prevbest, nvars)
double delta[VARS];
double point[VARS];
double prevbest;
{
    double minf;
    double z[VARS];
    double ftmp;

    int i;

    minf = prevbest;

    for (i = 0; i < nvars; i++)
    {
        z[i] = point[i];
    }

    for (i = 0; i < nvars; i++)
    {
        z[i] = point[i] + delta[i];

        ftmp = f(z, nvars);

        if (ftmp < minf)
        {
            minf = ftmp;
        }
        else
        {
            delta[i] = 0.0 - delta[i];
            z[i]     = point[i] + delta[i];

            ftmp = f(z, nvars);

            if (ftmp < minf)
            {
                minf = ftmp;
            }
            else
            {
                z[i] = point[i];
            }
        }
    }

    for (i = 0; i < nvars; i++)
    {
        point[i] = z[i];
    }

    return minf;
}

/**
 * Main optimization function.
 * <br />
 * <br />The hooke subroutine itself.
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
hooke(nvars, startpt, endpt, rho, epsilon, itermax)
double startpt[VARS];
double endpt[VARS];
double rho;
double epsilon;
{
    int i;
    int iadj;
    int iters;
    int j;
    int keep;

    double newx[VARS];
    double xbefore[VARS];
    double delta[VARS];
    double steplength;
    double fbefore;
    double newf;
    double tmp;

    for (i = 0; i < nvars; i++)
    {
        newx[i] = xbefore[i] = startpt[i];

        delta[i] = fabs(startpt[i] * rho);

        if (delta[i] == 0.0)
        {
            delta[i] = rho;
        }
    }

    iadj       = 0;
    steplength = rho;
    iters      = 0;

    fbefore = f(newx, nvars);

    newf = fbefore;

    while ((iters < itermax) && (steplength > epsilon))
    {
        iters++;
        iadj++;

        printf("\nAfter %5d funevals, f(x) =  %.4e at\n", FUNEVALS, fbefore);

        for (j = 0; j < nvars; j++)
        {
            printf("   x[%2d] = %.4e\n", j, xbefore[j]);
        }

        /* Find best new point, one coord at a time. */
        for (i = 0; i < nvars; i++)
        {
            newx[i] = xbefore[i];
        }

        newf = best_nearby(delta, newx, fbefore, nvars);

        /* If we made some improvements, pursue that direction. */
        keep = 1;

        while ((newf < fbefore) && (keep == 1))
        {
            iadj = 0;

            for (i = 0; i < nvars; i++)
            {
                /* Firstly, arrange the sign of delta[]. */
                if (newx[i] <= xbefore[i])
                {
                    delta[i] = 0.0 - fabs(delta[i]);
                }
                else
                {
                    delta[i] = fabs(delta[i]);
                }

                /* Now, move further in this direction. */
                tmp        = xbefore[i];
                xbefore[i] = newx[i];
                newx[i]    = newx[i] + newx[i] - tmp;
            }

            fbefore = newf;

            newf = best_nearby(delta, newx, fbefore, nvars);

            /* If the further (optimistic) move was bad.... */
            if (newf >= fbefore)
            {
                break;
            }

            /*
             * Make sure that the differences between the new and the old
             * points are due to actual displacements; beware of roundoff
             * errors that might cause newf < fbefore.
             */
            keep = 0;

            for (i = 0; i < nvars; i++)
            {
                keep = 1;

                if (fabs(newx[i] - xbefore[i]) > (0.5 * fabs(delta[i])))
                {
                    break;
                }
                else
                {
                    keep = 0;
                }
            }
        }

        if ((steplength >= epsilon) && (newf >= fbefore))
        {
            steplength = steplength * rho;

            for (i = 0; i < nvars; i++)
            {
                delta[i] *= rho;
            }
        }
    }

    for (i = 0; i < nvars; i++)
    {
        endpt[i] = xbefore[i];
    }

    return iters;
}

/**
 * The user-supplied objective function f(x,n).
 * <br />
 * <br />Represents either the Rosenbrock's classic
 * parabolic valley (&quot;banana&quot;) function
 * or the so-called &quot;Woods&quot; function.
 *
 * @param x The point at which f(x) should be evaluated.
 * @param n The number of coordinates of <code>x</code>.
 *
 * @return The objective function value.
 */
double
f(x, n)
double x[VARS];
{
#ifndef WOODS
    /* Rosenbrock's classic parabolic valley ("banana") function. */
    double a;
    double b;
    double c;

    FUNEVALS++;

    a = x[0];
    b = x[1];

    c = 100.0 * (b - (a * a)) * (b - (a * a));

    return (c + ((1.0 - a) * (1.0 - a)));
#else
    /* Woods -- a la More, Garbow & Hillstrom (TOMS algorithm 566). */
    double s1;
    double s2;
    double s3;
    double t1;
    double t2;
    double t3;
    double t4;
    double t5;

    FUNEVALS++;

    s1 = x[1] - x[0] * x[0];
    s2 = 1    - x[0];
    s3 = x[1] - 1;

    t1 = x[3] - x[2] * x[2];
    t2 = 1    - x[2];
    t3 = x[3] - 1;

    t4 = s3 + t3;
    t5 = s3 - t3;

    return (100 * (s1 * s1) + s2 * s2
           + 90 * (t1 * t1) + t2 * t2
           + 10 * (t4 * t4) + t5 * t5 / 10.);
#endif
}

/* Main program function main() :-). */
main()
{
    int nvars;
    int itermax;
    int jj;
    int i;

    double startpt[VARS];
    double rho;
    double epsilon;
    double endpt[VARS];

#ifndef WOODS
    /* Starting guess for Rosenbrock's test function. */
    nvars      =  2;
    startpt[0] = -1.2;
    startpt[1] =  1.0;
    rho        =  RHO_BEGIN;
#else
    /* Starting guess test problem "Woods". */
    nvars      =  4;
    startpt[0] = -3;
    startpt[1] = -1;
    startpt[2] = -3;
    startpt[3] = -1;
    rho        =  RHO_WOODS;
#endif

    itermax = IMAX;
    epsilon = EPSMIN;

    jj = hooke(nvars, startpt, endpt, rho, epsilon, itermax);

    printf("\n\n\nHOOKE USED %d ITERATIONS, AND RETURNED\n", jj);

    for (i = 0; i < nvars; i++)
    {
        printf("x[%3d] = %15.7e \n", i, endpt[i]);
    }

#ifdef WOODS
    puts("True answer: f(1, 1, 1, 1) = 0.");
#endif
}

/* vim:set nu:et:ts=4:sw=4: */
