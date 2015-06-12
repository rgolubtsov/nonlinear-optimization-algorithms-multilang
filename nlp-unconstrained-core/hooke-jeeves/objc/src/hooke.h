/*
 * nlp-unconstrained-core/hooke-jeeves/objc/src/hooke.h
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

#import <Foundation/Foundation.h>

/** Constant. The maximum number of variables. */
extern const NSUInteger VARS;

/** Constant. The ending value of stepsize. */
extern const NSUInteger EPSMIN;

/** Constant. The maximum number of iterations. */
extern const NSUInteger IMAX;

/** Helper constant. */
extern const NSUInteger INDEX_ZERO;

/** Helper constant. */
extern const NSUInteger INDEX_ONE;

/** Helper constant. */
extern const CGFloat ZERO_POINT_FIVE;

/**
 * The <code>Hooke</code> class contains methods for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 *
 * @author  Radislav (Radic) Golubtsov
 * @version 0.1
 * @see     Rosenbrock
 * @see     Woods
 * @since   hooke-jeeves 0.1
 */
@interface Hooke : NSObject
{
@private
    /*
     * GCC kludge: Properties must be declared as ivars too
     *             to avoid compile-time errors like the following:
     *
     * error: ivar ‘funEvals’ used by ‘@synthesize’ declaration
     *        must be an existing ivar
     *
     * Effective at least for GCC 5.1.0.
     */
    NSInteger funEvals;
}

/** The number of function evaluations. */
@property NSInteger funEvals;

/**
 * Helper method.
 * <br />
 * <br />Given a point, look for a better one nearby, one coord at a time.
 *
 * @param delta    The delta between <code>prevBest</code>
 *                 and <code>point</code>.
 * @param point    The coordinate from where to begin.
 * @param prevBest The previous best-valued coordinate.
 * @param nVars    The number of variables.
 *
 * @returns The objective function value at a nearby.
 */
- (CGFloat) bestNearby : (CGFloat *) delta
               point__ : (CGFloat *) point
            prevBest__ : (CGFloat) prevBest
               nVars__ : (NSUInteger) nVars;

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
 * @returns The number of iterations used to find the local minimum.
 */
- (NSUInteger) hooke : (NSUInteger) nVars
           startPt__ : (CGFloat *) startPt
             endPt__ : (CGFloat *) endPt
               rho__ : (CGFloat) rho
           epsilon__ : (CGFloat) epsilon
           iterMax__ : (NSUInteger) iterMax;

@end

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
