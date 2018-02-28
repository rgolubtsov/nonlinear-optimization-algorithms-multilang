/*
 * nlp-unconstrained-cli/hooke-jeeves/objc/src/rosenbrock.h
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

#import "hooke.h"

/** Helper constant. */
extern const CGFloat ONE_HUNDRED_POINT_ZERO;

/**
 * The <code>Rosenbrock</code> class is responsible for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 * <br />
 * <br />The objective function in this case
 * is the Rosenbrock's parabolic valley function.
 *
 * @author  Radislav (Radicchio) Golubtsov
 * @version 0.1
 * @see     Hooke
 * @since   hooke-jeeves 0.1
 */
@interface Rosenbrock : NSObject
/**
 * The user-supplied objective function f(x,n).
 * <br />
 * <br />Represents here the Rosenbrock's classic parabolic valley
 * (&quot;banana&quot;) function.
 *
 * @param x         The point at which f(x) should be evaluated.
 * @param n         The number of coordinates of <code>x</code>.
 * @param cFunEvals The number of function evaluations container (FunEvals *).
 *
 * @return The objective function value.
 */
+ (CGFloat) f : (CGFloat *) x
          n__ : (NSUInteger) n
  cFunEvals__ : (id) cFunEvals;

@end

// vim:set nu et ts=4 sw=4:
