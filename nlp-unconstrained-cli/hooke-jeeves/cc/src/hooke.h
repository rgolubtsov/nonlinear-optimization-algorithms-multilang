/*
 * nlp-unconstrained-cli/hooke-jeeves/cc/src/hooke.h
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * Written by Radislav (Radicchio) Golubtsov, 2016
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

#ifndef __CC__HOOKE_H
#define __CC__HOOKE_H

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

/**
 * The <code>NLPUCCLIHooke</code> namespace is used as a container
 * for the <code>Hooke</code> class.
 */
namespace NLPUCCLIHooke {

/** Constant. The maximum number of variables. */
extern const unsigned int VARS;

/** Constant. The stepsize geometric shrink. */
extern const double RHO_BEGIN;

/**
 * Constant. The stepsize geometric shrink.
 * <br />
 * <br />The Hooke &amp; Jeeves algorithm works reasonably well
 * on Rosenbrock's function, but can fare worse on some standard
 * test functions, depending on rho. Here is an example that works well
 * when rho = 0.5, but fares poorly with rho = 0.6, and better again
 * with rho = 0.8.
 */
extern const double RHO_WOODS;

/** Constant. The ending value of stepsize. */
extern const double EPSMIN;

/** Constant. The maximum number of iterations. */
extern const unsigned int IMAX;

/** Helper constants. */
extern const unsigned int INDEX_ZERO;
extern const unsigned int INDEX_ONE;
extern const unsigned int INDEX_TWO;
extern const unsigned int INDEX_THREE;
extern const unsigned int TWO;
extern const unsigned int FOUR;
extern const double       MINUS_ONE_POINT_TWO;
extern const double       ONE_POINT_ZERO;
extern const int          MINUS_THREE;
extern const int          MINUS_ONE;
extern const double       ZERO_POINT_FIVE;

/**
 * The <code>Hooke</code> class contains methods for solving a nonlinear
 * optimization problem using the algorithm of Hooke and Jeeves.
 *
 * @author  Radislav (Radicchio) Golubtsov
 * @version 0.1
 * @see     Rosenbrock
 * @see     Woods
 * @since   hooke-jeeves 0.1
 */
class Hooke {
private:
    /**
     * Helper method.
     * <br />
     * <br />Given a point, look for a better one nearby, one coord at a time.
     *
     * @param delta     The delta between <code>prevBest</code>
     *                  and <code>point</code>.
     * @param point     The coordinate from where to begin.
     * @param prevBest  The previous best-valued coordinate.
     * @param nVars     The number of variables.
     * @param cFunEvals The number of function evaluations container
     *                  (FunEvals *).
     *
     * @return The objective function value at a nearby.
     */
    double bestNearby(double *,
                      double *,
                      const double,
                      const unsigned int,
                      const void *);

public:
    /**
     * Main optimization method.
     * <br />
     * <br />The hooke subroutine itself.
     *
     * @param nVars   The number of variables.
     * @param startPt The starting point coordinates.
     * @param endPt   The ending point coordinates.
     * @param rho     The rho value.
     * @param epsilon The epsilon value.
     * @param iterMax The maximum number of iterations.
     *
     * @return The number of iterations used to find the local minimum.
     */
    unsigned int hooke(const unsigned int,
                       const double *,
                       double *,
                       const double,
                       const double,
                       const unsigned int);

    /** Default constructor. */
    Hooke();

    /** Destructor. */
    ~Hooke();
};

} // namespace NLPUCCLIHooke

#endif // __CC__HOOKE_H

// vim:set nu:et:ts=4:sw=4:
