/*
 * nlp-unconstrained-cli/hooke-jeeves/cc/src/rosenbrock.h
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

#ifndef __CC__ROSENBROCK_H
#define __CC__ROSENBROCK_H

#include "hooke.h"

/**
 * The <code>NLPUCCLIHooke</code> namespace is used as a container
 * for the <code>Rosenbrock</code> class.
 */
namespace NLPUCCLIHooke {

/** Helper constant. */
extern const double ONE_HUNDRED_POINT_ZERO;

/**
 * The <code>Rosenbrock</code> class is responsible for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 * <br />
 * <br />The objective function in this case
 * is the Rosenbrock's parabolic valley function.
 *
 * @author  Radislav (Radicchio) Golubtsov
 * @version 0.1
 * @see     Hooke
 * @since   hooke-jeeves 0.1
 */
class Rosenbrock {
public:
    /**
     * The user-supplied objective function f(x,n).
     * <br />
     * <br />Represents here the Rosenbrock's classic parabolic valley
     * (&quot;banana&quot;) function.
     *
     * @param x         The point at which f(x) should be evaluated.
     * @param n         The number of coordinates of <code>x</code>.
     * @param cFunEvals The number of function evaluations container
     *                  (FunEvals *).
     *
     * @return The objective function value.
     */
    static double f(const double *, const unsigned int, const void *);

    /** Default constructor. */
    Rosenbrock();

    /** Destructor. */
    ~Rosenbrock();
};

} // namespace NLPUCCLIHooke

#endif // __CC__ROSENBROCK_H

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
