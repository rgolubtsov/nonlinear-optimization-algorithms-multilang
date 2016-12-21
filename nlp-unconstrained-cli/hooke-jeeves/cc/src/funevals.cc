/*
 * nlp-unconstrained-core/hooke-jeeves/cc/src/funevals.cc
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

#include "funevals.h"

// The NLPUCCoreHooke namespace.
namespace NLPUCCoreHooke {

// Getter for funEvals.
unsigned int FunEvals::getFunEvals() {
    return funEvals;
}

// Setter for funEvals.
void FunEvals::setFunEvals(const unsigned int __funEvals) {
    funEvals = __funEvals;
}

// Default constructor.
FunEvals::FunEvals() {}

// Destructor.
FunEvals::~FunEvals() {}

} // namespace NLPUCCoreHooke

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
