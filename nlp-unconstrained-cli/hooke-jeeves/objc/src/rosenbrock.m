/*
 * nlp-unconstrained-cli/hooke-jeeves/objc/src/rosenbrock.m
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * Written by Radislav (Radicchio) Golubtsov, 2015-2020
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

#import "rosenbrock.h"
#import "funevals.h"

// Helper constant.
const CGFloat ONE_HUNDRED_POINT_ZERO = 100.0;

// The Rosenbrock class.
@implementation Rosenbrock
// The user-supplied objective function f(x,n).
+ (CGFloat) f : (CGFloat *) x
          n__ : (NSUInteger) n
  cFunEvals__ : (id) cFunEvals {

    CGFloat a;
    CGFloat b;
    CGFloat c;

    // Note: Since we're on Objective-C 2.0,
    //       using a "dot syntax" to access props.
//    [(FunEvals *) cFunEvals setFunEvals  :
//    [(FunEvals *) cFunEvals    funEvals] + 1];
    ((FunEvals *) cFunEvals).funEvals++;

    a = x[INDEX_ZERO];
    b = x[INDEX_ONE];

    c = ONE_HUNDRED_POINT_ZERO * (b - (a * a)) * (b - (a * a));

    return (c + ((ONE_POINT_ZERO - a) * (ONE_POINT_ZERO - a)));
}

@end

// vim:set nu et ts=4 sw=4:
