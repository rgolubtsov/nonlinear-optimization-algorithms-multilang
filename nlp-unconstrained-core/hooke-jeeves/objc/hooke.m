/*
 * nlp-unconstrained-core/hooke-jeeves/objc/hooke.m
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

#import "hooke.h"
#import "rosenbrock.h"
#import "woods.h"

#ifndef Woods

int main(void) {
    ORosenbrock *rosenbrock = [[ORosenbrock alloc] init];

    int i;

    // Starting guess for Rosenbrock's test function.
    rosenbrock->nvars      = 2;
    rosenbrock->startpt[0] = -1.2;
    rosenbrock->startpt[1] = 1.0;
    rosenbrock->itermax    = IMAX;
    rosenbrock->rho        = RHO_BEGIN;
    rosenbrock->epsilon    = EPSMIN;

    // Performing the Hooke-Jeeves search and printing the results.
    printf("\n\n\nHOOKE USED %d ITERATIONS, AND RETURNED\n",
           [rosenbrock hooke]);

    for (i = 0; i < rosenbrock->nvars; i++) {
        printf("x[%3d] = %15.7le \n", i, rosenbrock->endpt[i]);
    }

    return EXIT_SUCCESS;
}

#else

/**
 * The Hooke & Jeeves algorithm works reasonably well on Rosenbrock's
 * function, but can fare worse on some standard test functions,
 * depending on rho. Here is an example that works well when rho = 0.5,
 * but fares poorly with rho = 0.6, and better again with rho = 0.8.
 */
#ifndef RHO_WOODS
    #define RHO_WOODS 0.6
#endif

int main(void) {
    OWoods *woods = [[OWoods alloc] init];

    int i;

    // Starting guess test problem "Woods".
    woods->nvars      = 4;
    woods->startpt[0] = -3;
    woods->startpt[1] = -1;
    woods->startpt[2] = -3;
    woods->startpt[3] = -1;
    woods->itermax    = IMAX;
    woods->rho        = RHO_WOODS;
    woods->epsilon    = EPSMIN;

    // Performing the Hooke-Jeeves search and printing the results.
    printf("\n\n\nHOOKE USED %d ITERATIONS, AND RETURNED\n", [woods hooke]);

    for (i = 0; i < woods->nvars; i++) {
        printf("x[%3d] = %15.7le \n", i, woods->endpt[i]);
    }

    printf("True answer: f(1, 1, 1, 1) = 0.\n");

    return EXIT_SUCCESS;
}

#endif

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
