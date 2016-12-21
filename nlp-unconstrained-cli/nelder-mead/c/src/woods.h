/*
 * nlp-unconstrained-cli/nelder-mead/c/src/woods.h
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
 * The <code>woods.h</code> header file contains function prototypes
 * for solving a nonlinear optimization problem using the algorithm
 * of Nelder and Mead (&quot;Downhill simplex method&quot;).
 * <br />
 * <br />The objective function in this case
 * is the so-called &quot;Woods&quot; function.
 *
 * @author  Radislav (Radic) Golubtsov
 * @version 0.1
 * @since   nelder-mead 0.1
 */

#ifndef __C__WOODS_H
#define __C__WOODS_H

#include "nelmin.h"

/** Helper constants. */
#define ACOEFF 1
#define SCOEFF 100
#define TCOEFF 90
#define DCOEFF 10

/**
 * The user-supplied objective function f(x).
 * <br />
 * <br />Woods &ndash; a la More, Garbow &amp; Hillstrom
 * (TOMS algorithm 566).
 *
 * @param x The point at which f(x) should be evaluated.
 *
 * @return The objective function value.
 */
extern double f(const double *);

#endif /* __C__WOODS_H */

/* ========================================================================= */
/* vim:set nu:et:ts=4:sw=4:                                                  */
/* ========================================================================= */
