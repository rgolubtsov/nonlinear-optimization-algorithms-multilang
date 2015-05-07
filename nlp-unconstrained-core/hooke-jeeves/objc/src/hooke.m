/*
 * nlp-unconstrained-core/hooke-jeeves/objc/hooke.m
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * Copyright (C) 2015 Radislav (Radic) Golubtsov
 */

#import "hooke.h"
#import "rosenbrock.h"
#import "woods.h"

#ifndef Woods

int main(void) {
    ORosenbrock *rosenbrock = [[ORosenbrock alloc] init];

    int i;

    // Starting guess for Rosenbrock's test function.
    rosenbrock->nvars      = 2;
    rosenbrock->startpt[0] = -1.2;
    rosenbrock->startpt[1] = 1.0;
    rosenbrock->itermax    = IMAX;
    rosenbrock->rho        = RHO_BEGIN;
    rosenbrock->epsilon    = EPSMIN;

    // Performing the Hooke-Jeeves search and printing the results.
    printf("\n\n\nHOOKE USED %d ITERATIONS, AND RETURNED\n",
           [rosenbrock hooke]);

    for (i = 0; i < rosenbrock->nvars; i++) {
        printf("x[%3d] = %15.7le \n", i, rosenbrock->endpt[i]);
    }

    return EXIT_SUCCESS;
}

#else

/**
 * The Hooke & Jeeves algorithm works reasonably well on Rosenbrock's
 * function, but can fare worse on some standard test functions,
 * depending on rho. Here is an example that works well when rho = 0.5,
 * but fares poorly with rho = 0.6, and better again with rho = 0.8.
 */
#ifndef RHO_WOODS
    #define RHO_WOODS 0.6
#endif

int main(void) {
    OWoods *woods = [[OWoods alloc] init];

    int i;

    // Starting guess test problem "Woods".
    woods->nvars      = 4;
    woods->startpt[0] = -3;
    woods->startpt[1] = -1;
    woods->startpt[2] = -3;
    woods->startpt[3] = -1;
    woods->itermax    = IMAX;
    woods->rho        = RHO_WOODS;
    woods->epsilon    = EPSMIN;

    // Performing the Hooke-Jeeves search and printing the results.
    printf("\n\n\nHOOKE USED %d ITERATIONS, AND RETURNED\n", [woods hooke]);

    for (i = 0; i < woods->nvars; i++) {
        printf("x[%3d] = %15.7le \n", i, woods->endpt[i]);
    }

    printf("True answer: f(1, 1, 1, 1) = 0.\n");

    return EXIT_SUCCESS;
}

#endif

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
