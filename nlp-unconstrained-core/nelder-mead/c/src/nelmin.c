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
struct optimum *nelmin(const unsigned int  n,
                             double       *start,
                       const double        reqmin,
                       const double       *step,
                       const unsigned int  konvge,
                       const unsigned int  kcount) {

    struct optimum *opt; /* The structure to hold the optimum data          */
                         /* (and metadata) as the result of performing      */
                         /* the optimization procedure.                     */

    double xmin[VARS];   /* The coordinates of the point which is estimated */
                         /* to minimize the function.                       */

    double ynewlo;       /* The minimum value of the function.              */

    unsigned int icount; /* The number of function evaluations used.        */
    unsigned int numres; /* The number of restarts.                         */
    unsigned int ifault; /* The error indicator.                            */

    unsigned int jcount;
    unsigned int nn;
    unsigned int i;
    unsigned int j;
    unsigned int ilo;
    unsigned int ihi;
    unsigned int l;

    double dn;
    double dnn;
    double del;
    double rq;
    double p[VARS][VARS + 1];
    double y[VARS + 1];
    double x;
    double ylo;
    double z;
    double pbar[VARS];
    double pstar[VARS];
    double ystar;
    double p2star[VARS];
    double y2star;

    opt = malloc(sizeof(struct optimum));

    for (i = 0; i < n; i++) {
        opt->xmin[i] = 0;
    }

    opt->ynewlo = icount = numres = 0;

    opt->indics[INDEX_0] = icount;
    opt->indics[INDEX_1] = numres;

    /* Check the input parameters. */
    if (reqmin <= 0) {
        ifault = IFAULT_1;

        opt->indics[INDEX_2] = ifault;

        return (opt);
    }

    if (n < 1) {
        ifault = IFAULT_1;

        opt->indics[INDEX_2] = ifault;

        return (opt);
    }

    if (VARS < n) {
        ifault = IFAULT_1;

        opt->indics[INDEX_2] = ifault;

        return (opt);
    }

    if (konvge < 1) {
        ifault = IFAULT_1;

        opt->indics[INDEX_2] = ifault;

        return (opt);
    }

    jcount = konvge;
    dn     = n;
    nn     = n + 1;
    dnn    = nn;
    del    = 1;
    rq     = reqmin * dn;

    /* Construction of initial simplex. */
    L1000:;

    for (i = 0; i < n; i++) {
        p[i][nn] = start[i];
    }

    y[nn] = f(start);

    for (j = 0; j < n; j++) {
        x         = start[j];
        start[j] += step[j] * del;

        for (i = 0; i < n; i++) {
            p[i][j] = start[i];
        }

        y[j] = f(start);

        start[j] = x;
    }

    icount += nn;

    /* The simplex construction is complete. */

    /*
     * Find highest and lowest y values.
     * ynewlo = y[ihi] indicates the vertex of the simplex
     * to be replaced.
     */
    ylo = y[INDEX_0];
    ilo =   INDEX_0;

    for (i = 1; i < nn; i++) {
        if (y[i] < ylo) {
            ylo = y[i];
            ilo = i;
        }
    }

    L2000:;

    ynewlo = y[INDEX_0];
    ihi    =   INDEX_0;

    for (i = 1; i < nn; i++) {
        if (ynewlo < y[i]) {
            ynewlo = y[i];
            ihi    = i;
        }
    }

    /*
     * Calculate pbar, the centroid of the simplex vertices
     * excepting the vertex with y value ynewlo.
     */
    for (i = 0; i < n; i++) {
        z = 0;

        for (j = 0; j < nn; j++) {
            z += p[i][j];
        }

        z       -= p[i][ihi];
        pbar[i]  = z / dn;
    }

    /* Reflection through the centroid. */
    for (i = 0; i < n; i++) {
        pstar[i] = pbar[i] + RCOEFF * (pbar[i] - p[i][ihi]);
    }

    ystar = f(pstar);

    icount++;

    /* Successful reflection, so extension. */
    if (ystar < ylo) {
        for (i = 0; i < n; i++) {
            p2star[i] = pbar[i] + ECOEFF * (pstar[i] - pbar[i]);
        }

        y2star = f(p2star);

        icount++;

        /* Check extension. */
        if (ystar < y2star) {
            for (i = 0; i < n; i++) {
                p[i][ihi] = pstar[i];
            }

            y[ihi] = ystar;
        /* Retain extension or contraction. */
        } else {
            for (i = 0; i < n; i++) {
                p[i][ihi] = p2star[i];
            }

            y[ihi] = y2star;
        }
    /* No extension. */
    } else {
        l = 0;

        for (i = 0; i < nn; i++) {
            if (ystar < y[i]) {
                l++;
            }
        }

        if (1 < l) {
            for (i = 0; i < n; i++) {
                p[i][ihi] = pstar[i];
            }

            y[ihi] = ystar;
        /* Contraction on the y[ihi] side of the centroid. */
        } else if (l == 0) {
            for (i = 0; i < n; i++) {
                p2star[i] = pbar[i] + CCOEFF * (p[i][ihi] - pbar[i]);
            }

            y2star = f(p2star);

            icount++;

            /* Contract the whole simplex. */
            if (y[ihi] < y2star) {
                for (j = 0; j < nn; j++) {
                    for (i = 0; i < n; i++) {
                        p[i][j] = (p[i][j] + p[i][ilo]) * CCOEFF;
                        xmin[i] =  p[i][j];
                    }

                    y[j] = f(xmin);
                }

                icount += nn;

                if (kcount < icount) {
                    goto L3000;
                }

                ylo = y[INDEX_0];
                ilo =   INDEX_0;

                for (i = 1; i < nn; i++) {
                    if (y[i] < ylo) {
                        ylo = y[i];
                        ilo = i;
                    }
                }

                goto L2000;
            /* Retain contraction. */
            } else {
                for (i = 0; i < n; i++) {
                    p[i][ihi] = p2star[i];
                }

                y[ihi] = y2star;
            }
        /* Contraction on the reflection side of the centroid. */
        } else if (l == 1) {
            for (i = 0; i < n; i++) {
                p2star[i] = pbar[i] + CCOEFF * (pstar[i] - pbar[i]);
            }

            y2star = f(p2star);

            icount++;

            /* Retain reflection? */
            if (y2star <= ystar) {
                for (i = 0; i < n; i++) {
                    p[i][ihi] = p2star[i];
                }

                y[ihi] = y2star;
            } else {
                for (i = 0; i < n; i++) {
                    p[i][ihi] = pstar[i];
                }

                y[ihi] = ystar;
            }
        }
    }

    /* Check if ylo improved. */
    if (y[ihi] < ylo) {
        ylo = y[ihi];
        ilo = ihi;
    }

    jcount--;

    if (jcount != 0) {
        goto L2000;
    }

    /* Check to see if minimum reached. */
    if (icount <= kcount) {
        jcount = konvge;
        z      = 0;

        for (i = 0; i < nn; i++) {
            z += y[i];
        }

        x = z / dnn;
        z = 0;

        for (i = 0; i < nn; i++) {
            z += pow((y[i] - x), SQUARE);
        }

        if (rq < z) {
            goto L2000;
        }
    }

    /* Factorial tests to check that ynewlo is a local minimum. */
L3000:;

    for (i = 0; i < n; i++) {
        xmin[i] = p[i][ilo];
    }

    ynewlo = y[ilo];

    if (kcount < icount) {
        ifault = IFAULT_2;

        for (i = 0; i < n; i++) {
            opt->xmin[i] = xmin[i];
        }

        opt->ynewlo = ynewlo;

        opt->indics[INDEX_0] = icount;
        opt->indics[INDEX_1] = numres;
        opt->indics[INDEX_2] = ifault;

        return (opt);
    }

    ifault = IFAULT_0;

    for (i = 0; i < n; i++) {
        del      = step[i] * EPS;
        xmin[i] += del;

        z = f(xmin);

        icount++;

        if (z < ynewlo) {
            ifault = IFAULT_2;

            goto L4000;
        }

        xmin[i] -= del * 2;

        z = f(xmin);

        icount++;

        if (z < ynewlo) {
            ifault = IFAULT_2;

            goto L4000;
        }

        xmin[i] += del;
    }

L4000:;

    if (ifault == 0) {
        for (i = 0; i < n; i++) {
            opt->xmin[i] = xmin[i];
        }

        opt->ynewlo = ynewlo;

        opt->indics[INDEX_0] = icount;
        opt->indics[INDEX_1] = numres;
        opt->indics[INDEX_2] = ifault;

        return (opt);
    }

    /* Restart the procedure. */
    for (i = 0; i < n; i++) {
        start[i] = xmin[i];
    }

    del = EPS;

    numres++;

    goto L1000;
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

    struct optimum *opt;

#ifndef WOODS
    /* Starting guess for Rosenbrock's test function. */
    puts("\nTEST01\n  Apply NELMIN to ROSENBROCK function.");

    n              = ROSEN_GUESS_N;
    start[INDEX_0] = ROSEN_GUESS_1;
    start[INDEX_1] = ROSEN_GUESS_2;
#else
    /* Starting guess test problem "Woods". */
    puts("\nTEST05\n  Apply NELMIN to WOODS function.");

    n              = WOODS_GUESS_N;
    start[INDEX_0] = WOODS_GUESS_1;
    start[INDEX_1] = WOODS_GUESS_2;
    start[INDEX_2] = WOODS_GUESS_1;
    start[INDEX_3] = WOODS_GUESS_2;
#endif

    reqmin        = REQMIN_GUESS;
    step[INDEX_0] = STEP_GUESS_1;
    step[INDEX_1] = STEP_GUESS_2;
    konvge        = KONVGE_GUESS;
    kcount        = KCOUNT_GUESS;

    puts(  "\n  Starting point X:\n");

    for (i = 0; i < n; i++) {
        printf("  %14.6E\n", start[i]);
    }

    ynewlo = f(start);

    printf("\n  F(X) = %14.6E\n", ynewlo);

    opt = nelmin(n, start, reqmin, step, konvge, kcount);

    for (i = 0; i < n; i++) {
        xmin[i] = opt->xmin[i];
    }

    ynewlo = opt->ynewlo;

    icount = opt->indics[INDEX_0];
    numres = opt->indics[INDEX_1];
    ifault = opt->indics[INDEX_2];

    free(opt);

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
