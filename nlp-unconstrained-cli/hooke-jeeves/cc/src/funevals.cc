/*
 * nlp-unconstrained-cli/hooke-jeeves/cc/src/funevals.cc
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

#include "funevals.h"

// The NLPUCCLIHooke namespace.
namespace NLPUCCLIHooke {

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

} // namespace NLPUCCLIHooke

// vim:set nu et ts=4 sw=4:
