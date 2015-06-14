/*
 * nlp-unconstrained-core/hooke-jeeves/objc/src/funevals.h
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
 * The <code>FunEvals</code> class is a helper class.
 * It holds the only property &ndash; the number of objective function
 * evaluations (and corresponding generated accessor methods).
 *
 * @author  Radislav (Radic) Golubtsov
 * @version 0.1
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

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
