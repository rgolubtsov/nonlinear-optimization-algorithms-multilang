/*
 * nlp-unconstrained-cli/hooke-jeeves/c/src/hooke.c
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

#include "funevals.h"

#ifndef WOODS
    #include "rosenbrock.h"
#else
    #include "woods.h"
#endif

/* Helper function best_nearby(...). */
double best_nearby(double *delta,
                   double *point,
                   const double prevbest,
                   const unsigned int nvars,
                   void *__fun_evals) {

    double minf;
    double z[VARS];
    double ftmp;

    unsigned int i;

    minf = prevbest;

    for (i = 0; i < nvars; i++) {
        z[i] = point[i];
    }

    for (i = 0; i < nvars; i++) {
        z[i] = point[i] + delta[i];

#ifndef WOODS
        ftmp = f(z, nvars, __fun_evals);
#else
        ftmp = f(z, nvars, __fun_evals);
#endif

        if (ftmp < minf) {
            minf = ftmp;
        } else {
            delta[i] = 0.0 - delta[i];
            z[i]     = point[i] + delta[i];

#ifndef WOODS
            ftmp = f(z, nvars, __fun_evals);
#else
            ftmp = f(z, nvars, __fun_evals);
#endif

            if (ftmp < minf) {
                minf = ftmp;
            } else {
                z[i] = point[i];
            }
        }
    }

    for (i = 0; i < nvars; i++) {
        point[i] = z[i];
    }

    return minf;
}

/* Main optimization function hooke(...). */
unsigned int hooke(const unsigned int nvars,
                   const double *startpt,
                   double *endpt,
                   const double rho,
                   const double epsilon,
                   const unsigned int itermax) {

    unsigned int i;
    unsigned int iadj;
    unsigned int iters;
    unsigned int j;
    unsigned int keep;

    double newx[VARS];
    double xbefore[VARS];
    double delta[VARS];
    double steplength;
    double fbefore;
    double newf;
    double tmp;

    struct fun_evals *fe;

    for (i = 0; i < nvars; i++) {
        newx[i] = xbefore[i] = startpt[i];

        delta[i] = fabs(startpt[i] * rho);

        if (delta[i] == 0.0) {
            delta[i] = rho;
        }
    }

    iadj       = 0;
    steplength = rho;
    iters      = 0;

    /* Allocating memory for the fun_evals structure. */
    fe = malloc(sizeof(*fe));

#ifndef WOODS
    fbefore = f(newx, nvars, fe);
#else
    fbefore = f(newx, nvars, fe);
#endif

    newf = fbefore;

    while ((iters < itermax) && (steplength > epsilon)) {
        iters++;
        iadj++;

        printf(
            "\nAfter %5d funevals, f(x) =  %.4le at\n",
            get_funevals(fe), fbefore
        );

        for (j = 0; j < nvars; j++) {
            printf("   x[%2d] = %.4le\n", j, xbefore[j]);
        }

        /* Find best new point, one coord at a time. */
        for (i = 0; i < nvars; i++) {
            newx[i] = xbefore[i];
        }

        newf = best_nearby(delta, newx, fbefore, nvars, fe);

        /* If we made some improvements, pursue that direction. */
        keep = 1;

        while ((newf < fbefore) && (keep == 1)) {
            iadj = 0;

            for (i = 0; i < nvars; i++) {
                /* Firstly, arrange the sign of delta[]. */
                if (newx[i] <= xbefore[i]) {
                    delta[i] = 0.0 - fabs(delta[i]);
                } else {
                    delta[i] = fabs(delta[i]);
                }

                /* Now, move further in this direction. */
                tmp        = xbefore[i];
                xbefore[i] = newx[i];
                newx[i]    = newx[i] + newx[i] - tmp;
            }

            fbefore = newf;

            newf = best_nearby(delta, newx, fbefore, nvars, fe);

            /* If the further (optimistic) move was bad.... */
            if (newf >= fbefore) {
                break;
            }

            /*
             * Make sure that the differences between the new and the old
             * points are due to actual displacements; beware of roundoff
             * errors that might cause newf < fbefore.
             */
            keep = 0;

            for (i = 0; i < nvars; i++) {
                keep = 1;

                if (fabs(newx[i] - xbefore[i])
                    > (ZERO_POINT_FIVE * fabs(delta[i]))) {

                    break;
                } else {
                    keep = 0;
                }
            }
        }

        if ((steplength >= epsilon) && (newf >= fbefore)) {
            steplength = steplength * rho;

            for (i = 0; i < nvars; i++) {
                delta[i] *= rho;
            }
        }
    }

    for (i = 0; i < nvars; i++) {
        endpt[i] = xbefore[i];
    }

    /* Releasing memory, allocated for the fun_evals structure. */
    free(fe);

    return iters;
}

/* Main program function main() :-). */
int main(void) {
    unsigned int nvars;
    unsigned int itermax;
    unsigned int jj;
    unsigned int i;

    double startpt[VARS];
    double rho;
    double epsilon;
    double endpt[VARS];

#ifndef WOODS
    /* Starting guess for Rosenbrock's test function. */
    nvars                = TWO;
    startpt[INDEX_ZERO]  = MINUS_ONE_POINT_TWO;
    startpt[INDEX_ONE]   = ONE_POINT_ZERO;
    rho                  = RHO_BEGIN;
#else
    /* Starting guess test problem "Woods". */
    nvars                = FOUR;
    startpt[INDEX_ZERO]  = MINUS_THREE;
    startpt[INDEX_ONE]   = MINUS_ONE;
    startpt[INDEX_TWO]   = MINUS_THREE;
    startpt[INDEX_THREE] = MINUS_ONE;
    rho                  = RHO_WOODS;
#endif

    itermax = IMAX;
    epsilon = EPSMIN;

    jj = hooke(nvars, startpt, endpt, rho, epsilon, itermax);

    printf("\n\n\nHOOKE USED %d ITERATIONS, AND RETURNED\n", jj);

    for (i = 0; i < nvars; i++) {
        printf("x[%3d] = %15.7le \n", i, endpt[i]);
    }

#ifdef WOODS
    puts("True answer: f(1, 1, 1, 1) = 0.");
#endif

    return EXIT_SUCCESS;
}

/* vim:set nu:et:ts=4:sw=4: */
