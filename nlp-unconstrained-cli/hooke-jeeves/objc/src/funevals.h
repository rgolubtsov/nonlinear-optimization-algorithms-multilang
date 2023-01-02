/*
 * nlp-unconstrained-cli/hooke-jeeves/objc/src/funevals.h
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * Written by Radislav (Radicchio) Golubtsov, 2015-2023
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

/**
 * The <code>FunEvals</code> class is a helper class.
 * It holds the only property &ndash; the number of objective function
 * evaluations (and corresponding generated accessor methods).
 *
 * @author  Radislav (Radicchio) Golubtsov
 * @version 0.1.1
 * @see     Hooke
 * @since   hooke-jeeves 0.1
 */
@interface FunEvals : NSObject
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
    NSUInteger funEvals;
}

/** The number of function evaluations. */
@property NSUInteger funEvals;

@end

// vim:set nu et ts=4 sw=4:
