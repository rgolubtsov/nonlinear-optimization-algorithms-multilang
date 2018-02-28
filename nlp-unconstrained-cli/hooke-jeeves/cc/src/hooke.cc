/*
 * nlp-unconstrained-cli/hooke-jeeves/cc/src/hooke.cc
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * Written by Radislav (Radicchio) Golubtsov, 2015-2017
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

// The NLPUCCLIHooke namespace.
namespace NLPUCCLIHooke {

// Constant. The maximum number of variables.
const unsigned int VARS = 250;

// Constant. The stepsize geometric shrink.
const double RHO_BEGIN = 0.5;

// Constant. The stepsize geometric shrink.
const double RHO_WOODS = 0.6;

// Constant. The ending value of stepsize.
const double EPSMIN = 1E-6;

// Constant. The maximum number of iterations.
const unsigned int IMAX = 5000;

// Helper constants.
const unsigned int INDEX_ZERO          =  0;
const unsigned int INDEX_ONE           =  1;
const unsigned int INDEX_TWO           =  2;
const unsigned int INDEX_THREE         =  3;
const unsigned int TWO                 =  2;
const unsigned int FOUR                =  4;
const double       MINUS_ONE_POINT_TWO = -1.2;
const double       ONE_POINT_ZERO      =  1.0;
const int          MINUS_THREE         = -3;
const int          MINUS_ONE           = -1;
const double       ZERO_POINT_FIVE     =  0.5;

// Helper method bestNearby(...).
double Hooke::bestNearby(double *delta,
                         double *point,
                         const double prevBest,
                         const unsigned int nVars,
                         const void *cFunEvals) {

    double minF;
    double z[VARS];
    double fTmp;

    unsigned int i;

    minF = prevBest;

    for (i = 0; i < nVars; i++) {
        z[i] = point[i];
    }

    for (i = 0; i < nVars; i++) {
        z[i] = point[i] + delta[i];

#ifndef WOODS
        Rosenbrock r;

        fTmp = r.f(z, nVars, cFunEvals);
#else
        Woods w;

        fTmp = w.f(z, nVars, cFunEvals);
#endif

        if (fTmp < minF) {
            minF = fTmp;
        } else {
            delta[i] = 0.0 - delta[i];
            z[i]     = point[i] + delta[i];

#ifndef WOODS
            fTmp = r.f(z, nVars, cFunEvals);
#else
            fTmp = w.f(z, nVars, cFunEvals);
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
unsigned int Hooke::hooke(const unsigned int nVars,
                          const double *startPt,
                          double *endPt,
                          const double rho,
                          const double epsilon,
                          const unsigned int iterMax) {

    unsigned int i;
    unsigned int iAdj;
    unsigned int iters;
    unsigned int j;
    unsigned int keep;

    double newX[VARS];
    double xBefore[VARS];
    double delta[VARS];
    double stepLength;
    double fBefore;
    double newF;
    double tmp;

    for (i = 0; i < nVars; i++) {
        newX[i] = xBefore[i] = startPt[i];

        delta[i] = std::fabs(startPt[i] * rho);

        if (delta[i] == 0.0) {
            delta[i] = rho;
        }
    }

    iAdj       = 0;
    stepLength = rho;
    iters      = 0;

    // Instantiating the FunEvals class.
    FunEvals *fe = new FunEvals();

#ifndef WOODS
    Rosenbrock r;

    fBefore = r.f(newX, nVars, fe);
#else
    Woods w;

    fBefore = w.f(newX, nVars, fe);
#endif

    newF = fBefore;

    while ((iters < iterMax) && (stepLength > epsilon)) {
        iters++;
        iAdj++;

        std::cout << "\n" // Not using here std::endl -
                          // see http://en.cppreference.com/w/cpp/io/manip/endl
                          // for the reason why.
                  << "After " << std::setw(5) << fe->getFunEvals()
                  << " funevals, f(x) =  " << std::setprecision(4)
                  << std::scientific << fBefore << " at\n";

        for (j = 0; j < nVars; j++) {
            std::cout << "   x[" << std::setw(2) << j << "] = " << xBefore[j]
                      << "\n";
        }

        // Find best new point, one coord at a time.
        for (i = 0; i < nVars; i++) {
            newX[i] = xBefore[i];
        }

        newF = bestNearby(delta, newX, fBefore, nVars, fe);

        // If we made some improvements, pursue that direction.
        keep = 1;

        while ((newF < fBefore) && (keep == 1)) {
            iAdj = 0;

            for (i = 0; i < nVars; i++) {
                // Firstly, arrange the sign of delta[].
                if (newX[i] <= xBefore[i]) {
                    delta[i] = 0.0 - std::fabs(delta[i]);
                } else {
                    delta[i] = std::fabs(delta[i]);
                }

                // Now, move further in this direction.
                tmp        = xBefore[i];
                xBefore[i] = newX[i];
                newX[i]    = newX[i] + newX[i] - tmp;
            }

            fBefore = newF;

            newF = bestNearby(delta, newX, fBefore, nVars, fe);

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

                if (std::fabs(newX[i] - xBefore[i])
                    > (ZERO_POINT_FIVE * std::fabs(delta[i]))) {

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

    // Destroying the FunEvals class instance.
    delete fe;

    return iters;
}

// Default constructor.
Hooke::Hooke() {}

// Destructor.
Hooke::~Hooke() {}

} // namespace NLPUCCLIHooke

using namespace NLPUCCLIHooke;

// Main program function main() :-).
int main(void) {
    unsigned int nVars;
    unsigned int iterMax;
    unsigned int jj;
    unsigned int i;

    double startPt[VARS];
    double rho;
    double epsilon;
    double endPt[VARS];

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
    Hooke *h = new Hooke();

    jj = h->hooke(nVars, startPt, endPt, rho, epsilon, iterMax);

    std::cout << "\n\n\nHOOKE USED " << jj << " ITERATIONS, AND RETURNED\n";

    for (i = 0; i < nVars; i++) {
        std::cout << "x[" << std::setw(3) << i << "] = " << std::setw(15)
                  << std::setprecision(7) << endPt[i] << " \n";
    }

#ifdef WOODS
    std::cout << "True answer: f(1, 1, 1, 1) = 0." << std::endl;
#endif

    // Destroying the Hooke class instance.
    delete h;

    return EXIT_SUCCESS;
}

// vim:set nu et ts=4 sw=4:
