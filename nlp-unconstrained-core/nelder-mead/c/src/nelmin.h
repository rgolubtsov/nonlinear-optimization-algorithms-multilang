/*
 * nlp-unconstrained-core/nelder-mead/c/src/nelmin.h
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Nelder-Mead nonlinear unconstrained minimization algorithm.
 * ============================================================================
 */

/**
 * The <code>nelmin.h</code> header file contains function prototypes
 * for solving a nonlinear optimization problem using the algorithm
 * of Nelder and Mead (&quot;Downhill simplex method&quot;).
 *
 * @author  Radislav (Radic) Golubtsov
 * @version 0.1
 * @since   nelder-mead 0.1
 */

#ifndef __C__NELMIN_H
#define __C__NELMIN_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/** Constant. The maximum number of variables. */
#define VARS 20

/** Constant. The reflection coefficient. */
#define RCOEFF 1.0E+00

/** Constant. The extension coefficient. */
#define ECOEFF 2.0E+00

/** Constant. The contraction coefficient. */
#define CCOEFF 0.5E+00

/** Constant. The optimality factor. */
#define EPS 0.001E+00

/**
 * Main optimization function.
 * <br />
 * <br />The nelmin subroutine itself (Nelder-Mead minimization).
 *
 * @param f      The objective function f(x).
 * @param n      The number of variables.
 * @param start  The starting point for the iteration.
 * @param xmin   The coordinates of the point which is estimated
 *               to minimize the function.
 * @param ynewlo The minimum value of the function.
 * @param reqmin The terminating limit for the variance of function values.
 * @param step   The size and shape of the initial simplex.
 * @param konvge The convergence check.
 * @param kcount The maximum number of function evaluations.
 *
 * @return The array containing the following indicators:
 *         <ul><li>The number of function evaluations used.</li>
 *             <li>The number of restarts.</li>
 *             <li>The error indicator.</li></ul>
 */
extern unsigned int *nelmin(      double (*)(const double *),
                            const unsigned int,
                            const double *,
                            const double *,
                            const double,
                            const double,
                            const double *,
                            const          int,
                            const unsigned int);

#endif /* __C__NELMIN_H */

/* ========================================================================= */
/* vim:set nu:et:ts=4:sw=4:                                                  */
/* ========================================================================= */
