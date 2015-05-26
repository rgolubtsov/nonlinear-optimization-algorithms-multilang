/*
 * nlp-unconstrained-core/hooke-jeeves/objc/rosenbrock.m
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

#import "rosenbrock.h"

@implementation ORosenbrock
- (double) f: (double *) x atN: (int) n {
    double a;
    double b;
    double c;

    funevals++;

    a = x[0];
    b = x[1];

    c = 100.0 * (b - (a * a)) * (b - (a * a));

    return (c + ((1.0 - a) * (1.0 - a)));
}

- (double) bestNearby: (double *) delta
              atPoint: (double *) point
           atPrevbest: (double) prevbest {

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

        ftmp = [self f: z atN: nvars];

        if (ftmp < minf) {
            minf = ftmp;
        } else {
            delta[i] = 0.0 - delta[i];
            z[i]     = point[i] + delta[i];

            ftmp = [self f: z atN: nvars];

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

- (int) hooke {
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

    fbefore = [self f: newx atN: nvars];

    newf = fbefore;

    while ((iters < itermax) && (steplength > epsilon)) {
        iters++;
        iadj++;

        printf("\nAfter %5d funevals, f(x) =  %.4le at\n", funevals, fbefore);

        for (j = 0; j < nvars; j++) {
            printf("   x[%2d] = %.4le\n", j, xbefore[j]);
        }

        // Find best new point, one coord at a time.
        for (i = 0; i < nvars; i++) {
            newx[i] = xbefore[i];
        }

        newf = [self bestNearby: delta atPoint: newx atPrevbest: fbefore];

        // If we made some improvements, pursue that direction.
        keep = 1;

        while ((newf < fbefore) && (keep == 1)) {
            iadj = 0;

            for (i = 0; i < nvars; i++) {
                // Firstly, arrange the sign of delta[].
                if (newx[i] <= xbefore[i]) {
                    delta[i] = 0.0 - fabs(delta[i]);
                } else {
                    delta[i] = fabs(delta[i]);
                }

                // Now, move further in this direction.
                tmp        = xbefore[i];
                xbefore[i] = newx[i];
                newx[i]    = newx[i] + newx[i] - tmp;
            }

            fbefore = newf;

            newf = [self bestNearby: delta atPoint: newx atPrevbest: fbefore];

            // If the further (optimistic) move was bad....
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
@end

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
