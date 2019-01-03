/*
 * nlp-unconstrained-cli/hooke-jeeves/objc/src/hooke.h
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * Written by Radislav (Radicchio) Golubtsov, 2015-2019
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

// (Not so dirty :-) hack: But we need double! Even on 32-bit systems.
#define CGFLOAT_DEFINED 1
typedef double          CGFloat;
#define CGFLOAT_MIN     DBL_MIN
#define CGFLOAT_MAX     DBL_MAX
#define CGFLOAT_IS_DBL  1

#import <Foundation/Foundation.h>

/** Constant. The maximum number of variables. */
extern const NSUInteger VARS;

/** Constant. The stepsize geometric shrink. */
extern const CGFloat RHO_BEGIN;

/**
 * Constant. The stepsize geometric shrink.
 * <br />
 * <br />The Hooke &amp; Jeeves algorithm works reasonably well
 * on Rosenbrock's function, but can fare worse on some standard
 * test functions, depending on rho. Here is an example that works well
 * when rho = 0.5, but fares poorly with rho = 0.6, and better again
 * with rho = 0.8.
 */
extern const CGFloat RHO_WOODS;

/** Constant. The ending value of stepsize. */
extern const CGFloat EPSMIN;

/** Constant. The maximum number of iterations. */
extern const NSUInteger IMAX;

/** Helper constants. */
extern const NSUInteger INDEX_ZERO;
extern const NSUInteger INDEX_ONE;
extern const NSUInteger INDEX_TWO;
extern const NSUInteger INDEX_THREE;
extern const NSUInteger TWO;
extern const NSUInteger FOUR;
extern const CGFloat    MINUS_ONE_POINT_TWO;
extern const CGFloat    ONE_POINT_ZERO;
extern const NSInteger  MINUS_THREE;
extern const NSInteger  MINUS_ONE;
extern const CGFloat    ZERO_POINT_FIVE;

/**
 * The <code>Hooke</code> class contains methods for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 *
 * @author  Radislav (Radicchio) Golubtsov
 * @version 0.1
 * @see     Rosenbrock
 * @see     Woods
 * @since   hooke-jeeves 0.1
 */
@interface Hooke : NSObject
/**
 * Helper method.
 * <br />
 * <br />Given a point, look for a better one nearby, one coord at a time.
 *
 * @param delta     The delta between <code>prevBest</code>
 *                  and <code>point</code>.
 * @param point     The coordinate from where to begin.
 * @param prevBest  The previous best-valued coordinate.
 * @param nVars     The number of variables.
 * @param cFunEvals The number of function evaluations container (FunEvals *).
 *
 * @return The objective function value at a nearby.
 */
- (CGFloat) bestNearby : (CGFloat *) delta
               point__ : (CGFloat *) point
            prevBest__ : (CGFloat) prevBest
               nVars__ : (NSUInteger) nVars
           cFunEvals__ : (id) cFunEvals;

/**
 * Main optimization method.
 * <br />
 * <br />The hooke subroutine itself.
 *
 * @param nVars   The number of variables.
 * @param startPt The starting point coordinates.
 * @param endPt   The ending point coordinates.
 * @param rho     The rho value.
 * @param epsilon The epsilon value.
 * @param iterMax The maximum number of iterations.
 *
 * @return The number of iterations used to find the local minimum.
 */
- (NSUInteger) hooke : (NSUInteger) nVars
           startPt__ : (CGFloat *) startPt
             endPt__ : (CGFloat *) endPt
               rho__ : (CGFloat) rho
           epsilon__ : (CGFloat) epsilon
           iterMax__ : (NSUInteger) iterMax;

@end

// vim:set nu et ts=4 sw=4:
