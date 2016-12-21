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

/** Helper constants. */
extern const NSUInteger ONE_HUNDRED;
extern const NSUInteger NINETY;
extern const NSUInteger TEN;
extern const CGFloat    TEN_POINT;

/**
 * The <code>Woods</code> class is responsible for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 * <br />
 * <br />The objective function in this case
 * is the so-called &quot;Woods&quot; function.
 *
 * @author  Radislav (Radic) Golubtsov
 * @version 0.1
 * @see     Hooke
 * @since   hooke-jeeves 0.1
 */
@interface Woods : NSObject
/**
 * The user-supplied objective function f(x,n).
 * <br />
 * <br />Woods &ndash; a la More, Garbow &amp; Hillstrom
 * (TOMS algorithm 566).
 *
 * @param x         The point at which f(x) should be evaluated.
 * @param n         The number of coordinates of <code>x</code>.
 * @param cFunEvals The number of function evaluations container (FunEvals *).
 *
 * @return The objective function value.
 */
+ (CGFloat) f : (CGFloat *) x
          n__ : (NSUInteger) n
  cFunEvals__ : (id) cFunEvals;

@end

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
