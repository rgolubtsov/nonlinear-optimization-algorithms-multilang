/*
 * nlp-unconstrained-core/hooke-jeeves/cc/src/woods.h
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

#ifndef __CC__WOODS_H
#define __CC__WOODS_H

#include "hooke.h"

/**
 * The <code>NLPUCCoreHooke</code> namespace is used as a container
 * for the <code>Woods</code> class.
 */
namespace NLPUCCoreHooke {

/** Helper constant. */
extern const unsigned int ONE_HUNDRED;

/** Helper constant. */
extern const unsigned int NINETY;

/** Helper constant. */
extern const unsigned int TEN;

/** Helper constant. */
extern const double TEN_POINT;

/**
 * The <code>Woods</code> class is responsible for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 * <br />
 * <br />The objective function in this case
 * is so called &quot;Woods&quot; function.
 *
 * @author  Radislav (Radic) Golubtsov
 * @version 0.1
 * @see     Hooke
 * @since   hooke-jeeves 0.1
 */
class Woods {
public:
    /**
     * The user-supplied objective function f(x,n).
     * <br />
     * <br />Woods &ndash; a la More, Garbow &amp; Hillstrom
     * (TOMS algorithm 566).
     *
     * @param x       The point at which f(x) should be evaluated.
     * @param n       The number of coordinates of <code>x</code>.
     * @param fClsPtr The class instance containing the objective function.
     *
     * @returns The objective function value.
     */
    static double f(const double *, const unsigned int, const void *);

    /** Default constructor. */
    Woods();

    /** Destructor. */
    ~Woods();
};

} // namespace NLPUCCoreHooke

#endif // __CC__WOODS_H

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
