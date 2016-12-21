/*
 * nlp-unconstrained-core/hooke-jeeves/cc/src/rosenbrock.cc
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

#include "rosenbrock.h"
#include "funevals.h"

// The NLPUCCoreHooke namespace.
namespace NLPUCCoreHooke {

// Helper constant.
const double ONE_HUNDRED_POINT_ZERO = 100.0;

// The user-supplied objective function f(x,n).
double Rosenbrock::f(const double *x,
                     const unsigned int n,
                     const void *cFunEvals) {

    double a;
    double b;
    double c;

    ((FunEvals *) cFunEvals)->setFunEvals(
    ((FunEvals *) cFunEvals)->getFunEvals() + 1);

    a = x[INDEX_ZERO];
    b = x[INDEX_ONE];

    c = ONE_HUNDRED_POINT_ZERO * (b - (a * a)) * (b - (a * a));

    return (c + ((ONE_POINT_ZERO - a) * (ONE_POINT_ZERO - a)));
}

// Default constructor.
Rosenbrock::Rosenbrock() {}

// Destructor.
Rosenbrock::~Rosenbrock() {}

} // namespace NLPUCCoreHooke

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
