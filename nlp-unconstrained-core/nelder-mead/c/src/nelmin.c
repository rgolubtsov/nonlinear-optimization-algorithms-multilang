/*
 * nlp-unconstrained-core/nelder-mead/c/src/nelmin.c
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Nelder-Mead nonlinear unconstrained minimization algorithm.
 * ============================================================================
 */

#ifndef WOODS
    #include "rosenbrock.h"
#else
    #include "woods.h"
#endif

/* Main optimization function nelmin(...). */
unsigned int *nelmin(      double      (*f)(const double *),
                     const unsigned int  n,
                     const double       *start,
                     const double       *xmin,
                     const double        ynewlo,
                     const double        reqmin,
                     const double       *step,
                     const          int  konvge,
                     const unsigned int  kcount) {

    unsigned int *indics; /* The array containing the following indicators: */
    unsigned int  icount; /* The number of function evaluations used. */
    unsigned int  numres; /* The number of restarts. */
    unsigned int  ifault; /* The error indicator. */

    indics = malloc(sizeof(unsigned int) * 3);

    icount = numres = ifault = 0;

    indics[0] = icount;
    indics[1] = numres;
    indics[2] = ifault;

    /* TODO: Implement the subroutine body. */

    return (indics);
}

/* Main program function main() :-). */
int main(void) {
    unsigned int n;
    unsigned int konvge;
    unsigned int kcount;
    unsigned int i;
    unsigned int icount;
    unsigned int numres;
    unsigned int ifault;

    double start[VARS];
    double reqmin;
    double step[VARS];
    double ynewlo;
    double xmin[VARS];

    unsigned int *indics;

#ifndef WOODS
    /* Starting guess for Rosenbrock's test function. */
    puts("\nTEST01\n  Apply NELMIN to ROSENBROCK function.");

    n        =  2;
    start[0] = -1.2;
    start[1] =  1.0;
#else
    /* Starting guess test problem "Woods". */
    puts("\nTEST05\n  Apply NELMIN to WOODS function.");

    n        =  4;
    start[0] = -3.0;
    start[1] = -1.0;
    start[2] = -3.0;
    start[3] = -1.0;
#endif

    reqmin  = 1.0E-08;
    step[0] = 1.0;
    step[1] = 1.0;
    konvge  = 10;
    kcount  = 500;

    puts(  "\n  Starting point X:\n");

    for (i = 0; i < n; i++) {
        printf("  %14.6E\n", start[i]);
    }

    ynewlo = f(start);

    printf("\n  F(X) = %14.6E\n", ynewlo);

    indics = nelmin(f, n, start, xmin, ynewlo, reqmin, step, konvge, kcount);

    icount = indics[0];
    numres = indics[1];
    ifault = indics[2];

    free(indics);

    printf("\n  Return code IFAULT = %8i\n", ifault);

    puts(  "\n  Estimate of minimizing value X*:\n");

    for (i = 0; i < n; i++) {
        printf("  %14.6E\n", xmin[i]);
    }

    printf("\n  F(X*) = %14.6E\n", ynewlo);

    printf("\n  Number of iterations = %8i\n", icount);

    printf(  "  Number of restarts   = %8i\n", numres);

    return EXIT_SUCCESS;
}

/* ========================================================================= */
/* vim:set nu:et:ts=4:sw=4:                                                  */
/* ========================================================================= */
