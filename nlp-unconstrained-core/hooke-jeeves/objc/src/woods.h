/*
 * nlp-unconstrained-core/hooke-jeeves/objc/src/woods.h
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

#import "hooke.h"

/**
 * The Woods model.
 *
 * The Hooke-Jeeves nonlinear minimization algorithm in application
 * to the Woods test problem.
 */
@interface OWoods: NSObject {
@private
    int funevals;

@public
    int nvars;
    int itermax;

    double startpt[VARS];
    double rho;
    double epsilon;
    double endpt[VARS];
}

///---------------------
/// @name Test function
///---------------------

/**
 * The user-supplied objective function f(x,n).
 * Woods -- a la More, Garbow & Hillstrom (TOMS algorithm 566).
 *
 * @param x The point at which f(x) should be evaluated.
 * @param n The number of coordinates of x.
 *
 * @return The objective function value.
 */
- (double) f: (double *) x atN: (int) n;

///------------------------------
/// @name Hooke-Jeeves algorithm
///------------------------------

/**
 * Helper method.
 * Given a point, look for a better one nearby, one coord at a time.
 *
 * @param delta    The delta between prevbest and point.
 * @param point    The coordinate from where to begin.
 * @param prevbest The previous best-valued coordinate.
 *
 * @return The objective function value at a nearby.
 */
- (double) bestNearby: (double *) delta
              atPoint: (double *) point
           atPrevbest: (double) prevbest;

/**
 * Main optimization method. (See hooke.m for description.)
 *
 * @return The number of iterations used to find the local minimum.
 */
- (int) hooke;
@end

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
