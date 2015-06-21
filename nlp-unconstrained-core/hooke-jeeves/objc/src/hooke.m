/*
 * nlp-unconstrained-core/hooke-jeeves/objc/src/hooke.m
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

#import "hooke.h"
#import "funevals.h"

#ifndef WOODS
    #import "rosenbrock.h"
#else
    #import "woods.h"
#endif

// Constant. The maximum number of variables.
const NSUInteger VARS = 250;

// Constant. The stepsize geometric shrink.
const CGFloat RHO_BEGIN = 0.5;

// Constant. The stepsize geometric shrink.
const CGFloat RHO_WOODS = 0.6;

// Constant. The ending value of stepsize.
const CGFloat EPSMIN = 1E-6;

// Constant. The maximum number of iterations.
const NSUInteger IMAX = 5000;

// Helper constants.
const NSUInteger INDEX_ZERO          =  0;
const NSUInteger INDEX_ONE           =  1;
const NSUInteger INDEX_TWO           =  2;
const NSUInteger INDEX_THREE         =  3;
const NSUInteger TWO                 =  2;
const NSUInteger FOUR                =  4;
const CGFloat    MINUS_ONE_POINT_TWO = -1.2;
const CGFloat    ONE_POINT_ZERO      =  1.0;
const NSInteger  MINUS_THREE         = -3;
const NSInteger  MINUS_ONE           = -1;
const CGFloat    ZERO_POINT_FIVE     =  0.5;

// The Hooke class.
@implementation Hooke
// Helper method bestNearby(...).
- (CGFloat) bestNearby : (CGFloat *) delta
               point__ : (CGFloat *) point
            prevBest__ : (CGFloat) prevBest
               nVars__ : (NSUInteger) nVars
           cFunEvals__ : (id) cFunEvals {

    CGFloat minF;
    CGFloat z[VARS];
    CGFloat fTmp;

    NSUInteger i;

    minF = prevBest;

    for (i = 0; i < nVars; i++) {
        z[i] = point[i];
    }

    for (i = 0; i < nVars; i++) {
        z[i] = point[i] + delta[i];

#ifndef WOODS
        fTmp = [Rosenbrock f : z n__ : nVars cFunEvals__ : cFunEvals];
#else
        fTmp = [Woods      f : z n__ : nVars cFunEvals__ : cFunEvals];
#endif

        if (fTmp < minF) {
            minF = fTmp;
        } else {
            delta[i] = 0.0 - delta[i];
            z[i]     = point[i] + delta[i];

#ifndef WOODS
            fTmp = [Rosenbrock f : z n__ : nVars cFunEvals__ : cFunEvals];
#else
            fTmp = [Woods      f : z n__ : nVars cFunEvals__ : cFunEvals];
#endif

            if (fTmp < minF) {
                minF = fTmp;
            } else {
                z[i] = point[i];
            }
        }
    }

    for (i = 0; i < nVars; i++) {
        point[i] = z[i];
    }

    return minF;
}

// Main optimization method hooke(...).
- (NSUInteger) hooke : (NSUInteger) nVars
           startPt__ : (CGFloat *) startPt
             endPt__ : (CGFloat *) endPt
               rho__ : (CGFloat) rho
           epsilon__ : (CGFloat) epsilon
           iterMax__ : (NSUInteger) iterMax {

    NSUInteger i;
    NSUInteger iAdj;
    NSUInteger iters;
    NSUInteger j;
    NSUInteger keep;

    CGFloat newX[VARS];
    CGFloat xBefore[VARS];
    CGFloat delta[VARS];
    CGFloat stepLength;
    CGFloat fBefore;
    CGFloat newF;
    CGFloat tmp;

    for (i = 0; i < nVars; i++) {
        newX[i] = xBefore[i] = startPt[i];

        delta[i] = fabs(startPt[i] * rho);

        if (delta[i] == 0.0) {
            delta[i] = rho;
        }
    }

    iAdj       = 0;
    stepLength = rho;
    iters      = 0;

    // Instantiating the FunEvals class.
    // Note: Since there are no init args,
    //       (simplier) instantiating the class through "new".
//    FunEvals *fe = [[FunEvals alloc] init];
    FunEvals *fe = [FunEvals new];

#ifndef WOODS
    fBefore = [Rosenbrock f : newX n__ : nVars cFunEvals__ : fe];
#else
    fBefore = [Woods      f : newX n__ : nVars cFunEvals__ : fe];
#endif

    newF = fBefore;

    NSUInteger funEvals__;

    while ((iters < iterMax) && (stepLength > epsilon)) {
        iters++;
        iAdj++;

        // Note: Since we're on Objective-C 2.0,
        //       using a "dot syntax" to access props.
//        funEvals__ = [fe funEvals];
        funEvals__ = fe.funEvals;

        printf(
            "\nAfter %5d funevals, f(x) =  %.4le at\n",
            (unsigned int) funEvals__, fBefore
        );

        for (j = 0; j < nVars; j++) {
            printf("   x[%2d] = %.4le\n", (unsigned int) j, xBefore[j]);
        }

        // Find best new point, one coord at a time.
        for (i = 0; i < nVars; i++) {
            newX[i] = xBefore[i];
        }

        newF = [self bestNearby : delta
                        point__ : newX
                     prevBest__ : fBefore
                        nVars__ : nVars
                    cFunEvals__ : fe];

        // If we made some improvements, pursue that direction.
        keep = 1;

        while ((newF < fBefore) && (keep == 1)) {
            iAdj = 0;

            for (i = 0; i < nVars; i++) {
                // Firstly, arrange the sign of delta[].
                if (newX[i] <= xBefore[i]) {
                    delta[i] = 0.0 - fabs(delta[i]);
                } else {
                    delta[i] = fabs(delta[i]);
                }

                // Now, move further in this direction.
                tmp        = xBefore[i];
                xBefore[i] = newX[i];
                newX[i]    = newX[i] + newX[i] - tmp;
            }

            fBefore = newF;

            newF = [self bestNearby : delta
                            point__ : newX
                         prevBest__ : fBefore
                            nVars__ : nVars
                        cFunEvals__ : fe];

            // If the further (optimistic) move was bad....
            if (newF >= fBefore) {
                break;
            }

            /*
             * Make sure that the differences between the new and the old
             * points are due to actual displacements; beware of roundoff
             * errors that might cause newF < fBefore.
             */
            keep = 0;

            for (i = 0; i < nVars; i++) {
                keep = 1;

                if (fabs(newX[i] - xBefore[i])
                    > (ZERO_POINT_FIVE * fabs(delta[i]))) {

                    break;
                } else {
                    keep = 0;
                }
            }
        }

        if ((stepLength >= epsilon) && (newF >= fBefore)) {
            stepLength = stepLength * rho;

            for (i = 0; i < nVars; i++) {
                delta[i] *= rho;
            }
        }
    }

    for (i = 0; i < nVars; i++) {
        endPt[i] = xBefore[i];
    }

    return iters;
}

@end

// Main program function main() :-).
int main(void) {
    NSUInteger nVars;
    NSUInteger iterMax;
    NSUInteger jj;
    NSUInteger i;

    CGFloat startPt[VARS];
    CGFloat rho;
    CGFloat epsilon;
    CGFloat endPt[VARS];

#ifndef WOODS
    // Starting guess for Rosenbrock's test function.
    nVars                = TWO;
    startPt[INDEX_ZERO]  = MINUS_ONE_POINT_TWO;
    startPt[INDEX_ONE]   = ONE_POINT_ZERO;
    rho                  = RHO_BEGIN;
#else
    // Starting guess test problem "Woods".
    nVars                = FOUR;
    startPt[INDEX_ZERO]  = MINUS_THREE;
    startPt[INDEX_ONE]   = MINUS_ONE;
    startPt[INDEX_TWO]   = MINUS_THREE;
    startPt[INDEX_THREE] = MINUS_ONE;
    rho                  = RHO_WOODS;
#endif

    iterMax = IMAX;
    epsilon = EPSMIN;

    // Instantiating the Hooke class.
    // Note: Since there are no init args,
    //       (simplier) instantiating the class through "new".
//    Hooke *h = [[Hooke alloc] init];
    Hooke *h = [Hooke new];

    jj = [h hooke : nVars
        startPt__ : startPt
          endPt__ : endPt
            rho__ : rho
        epsilon__ : epsilon
        iterMax__ : iterMax];

    printf(
        "\n\n\nHOOKE USED %d ITERATIONS, AND RETURNED\n", (unsigned int) jj
    );

    for (i = 0; i < nVars; i++) {
        printf("x[%3d] = %15.7le \n", (unsigned int) i, endPt[i]);
    }

#ifdef WOODS
    printf("True answer: f(1, 1, 1, 1) = 0.\n");
#endif

    return EXIT_SUCCESS;
}

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
