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
struct optimum *nelmin(      double      (*f)(const double *),
                       const unsigned int  n,
                             double       *start,
                       const double        reqmin,
                       const double       *step,
                       const          int  konvge,
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

    opt->indics[0] = icount;
    opt->indics[1] = numres;

    /* Check the input parameters. */
    if (reqmin <= 0.0E+00) {
        ifault = 1;

        opt->indics[2] = ifault;

        return (opt);
    }

    if (n < 1) {
        ifault = 1;

        opt->indics[2] = ifault;

        return (opt);
    }

    if (VARS < n) {
        ifault = 1;

        opt->indics[2] = ifault;

        return (opt);
    }

    if (konvge < 1) {
        ifault = 1;

        opt->indics[2] = ifault;

        return (opt);
    }

    jcount = konvge;
    dn     = n;
    nn     = n + 1;
    dnn    = nn;
    del    = 1.0E+00;
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
    ylo = y[0];
    ilo = 1;

    for (i = 1; i < nn; i++) {
        if (y[i] < ylo) {
            ylo = y[i];
            ilo = i;
        }
    }

    L2000:;

    ynewlo = y[0];
    ihi    = 1;

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
        z = 0.0E+00;

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
                        p[i][j] = (p[i][j] + p[i][ilo]) * 0.5E+00;
                        xmin[i] =  p[i][j];
                    }

                    y[j] = f(xmin);
                }

                icount += nn;

                if (kcount < icount) {
                    goto L3000;
                }

                ylo = y[0];
                ilo = 1;

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
        z      = 0.0E+00;

        for (i = 0; i < nn; i++) {
            z += y[i];
        }

        x = z / dnn;
        z = 0.0E+00;

        for (i = 0; i < nn; i++) {
            z += pow((y[i] - x), 2);
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
        ifault = 2;

        for (i = 0; i < n; i++) {
            opt->xmin[i] = xmin[i];
        }

        opt->ynewlo = ynewlo;

        opt->indics[0] = icount;
        opt->indics[1] = numres;
        opt->indics[2] = ifault;

        return (opt);
    }

    ifault = 0;

    for (i = 0; i < n; i++) {
        del     = step[i] * EPS;
        xmin[i] += del;

        z = f(xmin);

        icount++;

        if (z < ynewlo) {
            ifault = 2;

            goto L4000;
        }

        xmin[i] -= del * 2;

        z = f(xmin);

        icount++;

        if (z < ynewlo) {
            ifault = 2;

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

        opt->indics[0] = icount;
        opt->indics[1] = numres;
        opt->indics[2] = ifault;

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

    opt = nelmin(f, n, start, reqmin, step, konvge, kcount);

    for (i = 0; i < n; i++) {
        xmin[i] = opt->xmin[i];
    }

    ynewlo = opt->ynewlo;

    icount = opt->indics[0];
    numres = opt->indics[1];
    ifault = opt->indics[2];

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
