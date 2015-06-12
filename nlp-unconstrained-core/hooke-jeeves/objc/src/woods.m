/*
 * nlp-unconstrained-core/hooke-jeeves/objc/src/woods.m
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

#import "woods.h"

// Constant. The stepsize geometric shrink.
const CGFloat RHO_WOODS = 0.6;

// Helper constant.
const NSUInteger INDEX_TWO = 2;

// Helper constant.
const NSUInteger INDEX_THREE = 3;

// Helper constant.
const NSUInteger ONE_HUNDRED = 100;

// Helper constant.
const NSUInteger NINETY = 90;

// Helper constant.
const NSUInteger TEN = 10;

// Helper constant.
const CGFloat TEN_POINT = 10.;

// Helper constant.
const NSUInteger FOUR = 4;

// Helper constant.
const NSInteger MINUS_THREE = -3;

// Helper constant.
const NSInteger MINUS_ONE = -1;

// The Woods class.
@implementation Woods
// The user-supplied objective function f(x,n).
+ (CGFloat) f : (CGFloat *) x n__ : (NSUInteger) n {
    CGFloat s1;
    CGFloat s2;
    CGFloat s3;
    CGFloat t1;
    CGFloat t2;
    CGFloat t3;
    CGFloat t4;
    CGFloat t5;

    // Instantiating the Hooke class.
    Hooke *h = [[Hooke alloc] init];

    [h setFunEvals : [h funEvals] + 1];

    s1 = x[INDEX_ONE] - x[INDEX_ZERO] * x[INDEX_ZERO];
    s2 = 1 - x[INDEX_ZERO];
    s3 = x[INDEX_ONE] - 1;
    t1 = x[INDEX_THREE] - x[INDEX_TWO] * x[INDEX_TWO];
    t2 = 1 - x[INDEX_TWO];
    t3 = x[INDEX_THREE] - 1;
    t4 = s3 + t3;
    t5 = s3 - t3;

    return (ONE_HUNDRED * (s1 * s1) + s2 * s2
               + NINETY * (t1 * t1) + t2 * t2
                  + TEN * (t4 * t4) + t5 * t5 / TEN_POINT);
}

@end

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
