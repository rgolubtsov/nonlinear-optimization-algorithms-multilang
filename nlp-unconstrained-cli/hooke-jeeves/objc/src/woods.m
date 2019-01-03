/*
 * nlp-unconstrained-cli/hooke-jeeves/objc/src/woods.m
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

#import "woods.h"
#import "funevals.h"

// Helper constants.
const NSUInteger ONE_HUNDRED = 100;
const NSUInteger NINETY      = 90;
const NSUInteger TEN         = 10;
const CGFloat    TEN_POINT   = 10.;

// The Woods class.
@implementation Woods
// The user-supplied objective function f(x,n).
+ (CGFloat) f : (CGFloat *) x
          n__ : (NSUInteger) n
  cFunEvals__ : (id) cFunEvals {

    CGFloat s1;
    CGFloat s2;
    CGFloat s3;
    CGFloat t1;
    CGFloat t2;
    CGFloat t3;
    CGFloat t4;
    CGFloat t5;

    // Note: Since we're on Objective-C 2.0,
    //       using a "dot syntax" to access props.
//    [(FunEvals *) cFunEvals setFunEvals  :
//    [(FunEvals *) cFunEvals    funEvals] + 1];
    ((FunEvals *) cFunEvals).funEvals++;

    s1 = x[INDEX_ONE]   - x[INDEX_ZERO] * x[INDEX_ZERO];
    s2 = 1              - x[INDEX_ZERO];
    s3 = x[INDEX_ONE]   - 1;

    t1 = x[INDEX_THREE] - x[INDEX_TWO]  * x[INDEX_TWO];
    t2 = 1              - x[INDEX_TWO];
    t3 = x[INDEX_THREE] - 1;

    t4 = s3 + t3;
    t5 = s3 - t3;

    return (ONE_HUNDRED * (s1 * s1) + s2 * s2
               + NINETY * (t1 * t1) + t2 * t2
                  + TEN * (t4 * t4) + t5 * t5 / TEN_POINT);
}

@end

// vim:set nu et ts=4 sw=4:
