/*
 * nlp-unconstrained-cli/nelder-mead/c/src/rosenbrock.h
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Nelder-Mead nonlinear unconstrained minimization algorithm.
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
 * The <code>rosenbrock.h</code> header file contains function prototypes
 * for solving a nonlinear optimization problem using the algorithm
 * of Nelder and Mead (&quot;Downhill simplex method&quot;).
 * <br />
 * <br />The objective function in this case
 * is the Rosenbrock's parabolic valley function.
 *
 * @author  Radislav (Radicchio) Golubtsov
 * @version 0.1
 * @since   nelder-mead 0.1
 */

#ifndef __C__ROSENBROCK_H
#define __C__ROSENBROCK_H

#include "nelmin.h"

/** Helper constants. */
#define ACOEFF 1
#define BCOEFF 100

/**
 * The user-supplied objective function f(x).
 * <br />
 * <br />Represents here the Rosenbrock's classic parabolic valley
 * (&quot;banana&quot;) function.
 *
 * @param x The point at which f(x) should be evaluated.
 *
 * @return The objective function value.
 */
extern double f(const double *);

#endif /* __C__ROSENBROCK_H */

/* vim:set nu et ts=4 sw=4: */
