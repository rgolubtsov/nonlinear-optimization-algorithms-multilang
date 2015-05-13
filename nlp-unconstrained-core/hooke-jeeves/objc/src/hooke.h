/*
 * nlp-unconstrained-core/hooke-jeeves/objc/hooke.h
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * Written by Radislav (Radic) Golubtsov <https://github.com/rgolubtsov>
 */

#import <Foundation/Foundation.h>
#import <math.h>

/// Max # of variables.
#define VARS 250

/// Stepsize geometric shrink.
#define RHO_BEGIN 0.5

/// Ending value of stepsize.
#define EPSMIN 1E-6

/// Max # of iterations.
#define IMAX 5000

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
