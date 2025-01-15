/*
 * nlp-unconstrained-cli/hooke-jeeves/cc/src/woods.h
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * Written by Radislav (Radicchio) Golubtsov, 2015-2025
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

#ifndef __CC__WOODS_H
#define __CC__WOODS_H

#include "hooke.h"

/**
 * The <code>NLPUCCLIHooke</code> namespace is used as a container
 * for the <code>Woods</code> class.
 */
namespace NLPUCCLIHooke {

/** Helper constants. */
extern const unsigned int ONE_HUNDRED;
extern const unsigned int NINETY;
extern const unsigned int TEN;
extern const double       TEN_POINT;

/**
 * The <code>Woods</code> class is responsible for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 * <br />
 * <br />The objective function in this case
 * is the so-called &quot;Woods&quot; function.
 *
 * @author  Radislav (Radicchio) Golubtsov
 * @version 0.1.1
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
     * @param x         The point at which f(x) should be evaluated.
     * @param n         The number of coordinates of <code>x</code>.
     * @param cFunEvals The number of function evaluations container
     *                  (FunEvals *).
     *
     * @return The objective function value.
     */
    static double f(const double *, const unsigned int, const void *);

    /** Default constructor. */
    Woods();

    /** Destructor. */
    ~Woods();
};

} // namespace NLPUCCLIHooke

#endif // __CC__WOODS_H

// vim:set nu et ts=4 sw=4:
