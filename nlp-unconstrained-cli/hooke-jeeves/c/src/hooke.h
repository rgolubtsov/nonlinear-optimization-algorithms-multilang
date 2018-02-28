/*
 * nlp-unconstrained-cli/hooke-jeeves/c/src/hooke.h
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * Written by Radislav (Radicchio) Golubtsov, 2015-2017
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

/**
 * The <code>hooke.h</code> header file contains function prototypes
 * for solving a nonlinear optimization problem using the algorithm
 * of Hooke and Jeeves.
 *
 * @author  Radislav (Radicchio) Golubtsov
 * @version 0.1
 * @since   hooke-jeeves 0.1
 */

#ifndef __C__HOOKE_H
#define __C__HOOKE_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/** Constant. The maximum number of variables. */
#define VARS 250

/** Constant. The stepsize geometric shrink. */
#define RHO_BEGIN 0.5

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

/** Constant. The ending value of stepsize. */
#define EPSMIN 1E-6

/** Constant. The maximum number of iterations. */
#define IMAX 5000

/** Helper constants. */
#define INDEX_ZERO           0
#define INDEX_ONE            1
#define INDEX_TWO            2
#define INDEX_THREE          3
#define TWO                  2
#define FOUR                 4
#define MINUS_ONE_POINT_TWO -1.2
#define ONE_POINT_ZERO       1.0
#define MINUS_THREE         -3
#define MINUS_ONE           -1
#define ZERO_POINT_FIVE      0.5

/**
 * Helper function.
 * <br />
 * <br />Given a point, look for a better one nearby, one coord at a time.
 *
 * @param delta       The delta between <code>prevbest</code>
 *                    and <code>point</code>.
 * @param point       The coordinate from where to begin.
 * @param prevbest    The previous best-valued coordinate.
 * @param nvars       The number of variables.
 * @param __fun_evals The number of function evaluations container (struct *).
 *
 * @return The objective function value at a nearby.
 */
extern double best_nearby(double *,
                          double *,
                          const double,
                          const unsigned int,
                          void *);

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
extern unsigned int hooke(const unsigned int,
                          const double *,
                          double *,
                          const double,
                          const double,
                          const unsigned int);

#endif /* __C__HOOKE_H */

/* vim:set nu et ts=4 sw=4: */
