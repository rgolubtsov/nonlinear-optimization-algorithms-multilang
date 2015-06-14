/*
 * nlp-unconstrained-core/hooke-jeeves/objc/src/rosenbrock.h
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

/** Helper constant. */
extern const CGFloat ONE_HUNDRED_POINT_ZERO;

/**
 * The <code>Rosenbrock</code> class is responsible for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 * <br />
 * <br />The objective function in this case
 * is the Rosenbrock's parabolic valley function.
 *
 * @author  Radislav (Radic) Golubtsov
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
 * @param x       The point at which f(x) should be evaluated.
 * @param n       The number of coordinates of <code>x</code>.
 * @param fClsPtr The class instance containing the objective function.
 *
 * @returns The objective function value.
 */
+ (CGFloat) f : (CGFloat *) x n__ : (NSUInteger) n fClsPtr__ : (id) fClsPtr;

@end

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
