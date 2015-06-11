/*
 * nlp-unconstrained-core/hooke-jeeves/c/src/hooke.c
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Max # of variables. */
#define VARS 250

/* Stepsize geometric shrink. */
#define RHO_BEGIN 0.5

/* Ending value of stepsize. */
#define EPSMIN 1E-6

/* Max # of iterations. */
#define IMAX 5000

/* Global variables. */
int funevals = 0;

#ifdef WOODS
    double f(double *, int);
#else
    /* Rosenbrock's classic parabolic valley ("banana") function. */
    double f(double *x, int n) {
        double a;
        double b;
        double c;

        funevals++;

        a = x[0];
        b = x[1];

        c = 100.0 * (b - (a * a)) * (b - (a * a));

        return (c + ((1.0 - a) * (1.0 - a)));
    }
#endif

/* Given a point, look for a better one nearby, one coord at a time. */
double best_nearby(double *delta, double *point, double prevbest, int nvars) {
    double minf;
    double z[VARS];
    double ftmp;

    int i;

    minf = prevbest;

    for (i = 0; i < nvars; i++) {
        z[i] = point[i];
    }

    for (i = 0; i < nvars; i++) {
        z[i] = point[i] + delta[i];

        ftmp = f(z, nvars);

        if (ftmp < minf) {
            minf = ftmp;
        } else {
            delta[i] = 0.0 - delta[i];
            z[i]     = point[i] + delta[i];

            ftmp = f(z, nvars);

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

int hooke(int nvars,
          double *startpt,
          double *endpt,
          double rho,
          double epsilon,
          int itermax) {

    int i;
    int iadj;
    int iters;
    int j;
    int keep;

    double newx[VARS];
    double xbefore[VARS];
    double delta[VARS];
    double steplength;
    double fbefore;
    double newf;
    double tmp;

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

    fbefore = f(newx, nvars);

    newf = fbefore;

    while ((iters < itermax) && (steplength > epsilon)) {
        iters++;
        iadj++;

        printf("\nAfter %5d funevals, f(x) =  %.4le at\n", funevals, fbefore);

        for (j = 0; j < nvars; j++) {
            printf("   x[%2d] = %.4le\n", j, xbefore[j]);
        }

        /* Find best new point, one coord at a time. */
        for (i = 0; i < nvars; i++) {
            newx[i] = xbefore[i];
        }

        newf = best_nearby(delta, newx, fbefore, nvars);

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

            newf = best_nearby(delta, newx, fbefore, nvars);

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

                if (fabs(newx[i] - xbefore[i]) > (0.5 * fabs(delta[i]))) {
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

    return iters;
}

#ifndef WOODS
    int main(void) {
        int nvars;
        int itermax;
        int jj;
        int i;

        double startpt[VARS];
        double rho;
        double epsilon;
        double endpt[VARS];

        /* Starting guess for Rosenbrock's test function. */
        nvars      = 2;
        startpt[0] = -1.2;
        startpt[1] = 1.0;
        itermax    = IMAX;
        rho        = RHO_BEGIN;
        epsilon    = EPSMIN;

        jj = hooke(nvars, startpt, endpt, rho, epsilon, itermax);

        printf("\n\n\nHOOKE USED %d ITERATIONS, AND RETURNED\n", jj);

        for (i = 0; i < nvars; i++) {
            printf("x[%3d] = %15.7le \n", i, endpt[i]);
        }

        return EXIT_SUCCESS;
    }
#else
    /*
     * The Hooke & Jeeves algorithm works reasonably well on Rosenbrock's
     * function, but can fare worse on some standard test functions,
     * depending on rho. Here is an example that works well when rho = 0.5,
     * but fares poorly with rho = 0.6, and better again with rho = 0.8.
     */
    #ifndef RHO_WOODS
        #define RHO_WOODS 0.6
    #endif

    /* Woods -- a la More, Garbow & Hillstrom (TOMS algorithm 566). */
    double f(double *x, int n) {
        double s1;
        double s2;
        double s3;
        double t1;
        double t2;
        double t3;
        double t4;
        double t5;

        funevals++;

        s1 = x[1] - x[0] * x[0];
        s2 = 1 - x[0];
        s3 = x[1] - 1;
        t1 = x[3] - x[2] * x[2];
        t2 = 1 - x[2];
        t3 = x[3] - 1;
        t4 = s3 + t3;
        t5 = s3 - t3;

        return (100 * (s1 * s1) + s2 * s2
               + 90 * (t1 * t1) + t2 * t2
               + 10 * (t4 * t4) + t5 * t5 / 10.);
    }

    int main(void) {
        int nvars;
        int itermax;
        int jj;
        int i;

        double startpt[VARS];
        double rho;
        double epsilon;
        double endpt[VARS];

        /* Starting guess test problem "Woods". */
        nvars      = 4;
        startpt[0] = -3;
        startpt[1] = -1;
        startpt[2] = -3;
        startpt[3] = -1;
        itermax    = IMAX;
        rho        = RHO_WOODS;
        epsilon    = EPSMIN;

        jj = hooke(nvars, startpt, endpt, rho, epsilon, itermax);

        printf("\n\n\nHOOKE USED %d ITERATIONS, AND RETURNED\n", jj);

        for (i = 0; i < nvars; i++) {
            printf("x[%3d] = %15.7le \n", i, endpt[i]);
        }

        printf("True answer: f(1, 1, 1, 1) = 0.\n");

        return EXIT_SUCCESS;
    }
#endif

/* ========================================================================= */
/* vim:set nu:et:ts=4:sw=4:                                                  */
/* ========================================================================= */
