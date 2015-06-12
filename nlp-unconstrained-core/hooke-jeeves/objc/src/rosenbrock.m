/*
 * nlp-unconstrained-core/hooke-jeeves/objc/src/rosenbrock.m
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

#import "rosenbrock.h"

// Constant. The stepsize geometric shrink.
const CGFloat RHO_BEGIN = 0.5;

// Helper constant.
const CGFloat ONE_HUNDRED_POINT_ZERO = 100.0;

// Helper constant.
const CGFloat ONE_POINT_ZERO = 1.0;

// Helper constant.
const NSUInteger TWO = 2;

// Helper constant.
const CGFloat MINUS_ONE_POINT_TWO = -1.2;

// The Rosenbrock class.
@implementation Rosenbrock
// The user-supplied objective function f(x,n).
+ (CGFloat) f : (CGFloat *) x n__ : (NSUInteger) n {
    CGFloat a;
    CGFloat b;
    CGFloat c;

    // Instantiating the Hooke class.
    Hooke *h = [[Hooke alloc] init];

    [h setFunEvals : [h funEvals] + 1];

    a = x[INDEX_ZERO];
    b = x[INDEX_ONE];

    c = ONE_HUNDRED_POINT_ZERO * (b - (a * a)) * (b - (a * a));

    return (c + ((ONE_POINT_ZERO - a) * (ONE_POINT_ZERO - a)));
}

@end

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
