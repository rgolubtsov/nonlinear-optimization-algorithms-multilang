/*
 * nlp-unconstrained-cli/nelder-mead/__orig/src/nelmin.txt.c
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Nelder-Mead nonlinear unconstrained minimization algorithm.
 * ============================================================================
 * (Re-)Written by Radislav (Radicchio) Golubtsov, 2015-2020
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

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * NELMIN minimizes a function using the Nelder-Mead algorithm.              *
 *                                                                           *
 * Discussion:                                                               *
 *   This routine seeks the minimum value of a user-specified function.      *
 *                                                                           *
 *   Simplex function minimization  procedure  due  to  Nelder+Mead  (1965), *
 *   as  implemented  by  O'Neill  (1971,  Appl.   Statist.   20,   338-45), *
 *   with  subsequent  comments  by  Chambers+Ertel   (1974,   23,   250-1), *
 *   Benyon (1976, 25, 97) and Hill (1978, 27, 380-2).                       *
 *                                                                           *
 *   The function to be minimized must be defined by a function of the form: *
 *                                                                           *
 *     FUNCTION F(X)                                                         *
 *     DOUBLE PRECISION F                                                    *
 *     DOUBLE PRECISION X(*)                                                 *
 *                                                                           *
 *   and  the  name  of  this   subroutine   must   be   declared   EXTERNAL *
 *   in the calling routine and passed as the argument F.                    *
 *                                                                           *
 *   This routine does not include a  termination  test  using  the  fitting *
 *   of a quadratic surface.                                                 *
 *                                                                           *
 * Modified:                                                                 *
 *   27 February 2008                                                        *
 *                                                                           *
 * Author:                                                                   *
 *   FORTRAN 77 version by R. O'Neill                                        *
 *   Modifications by John Burkardt                                          *
 *                                                                           *
 * Reference:                                                                *
 *   John Nelder, Roger Mead,                                                *
 *   A simplex method for function minimization,                             *
 *   Computer Journal, Volume 7, 1965, pages 308-313.                        *
 *                                                                           *
 *   R. O'Neill,                                                             *
 *   Algorithm AS 47:                                                        *
 *   Function Minimization Using a Simplex Procedure,                        *
 *   Applied Statistics, Volume 20, Number 3, 1971, pages 338-345.           *
 *                                                                           *
 * Parameters:                                                               *
 *   Input,  EXTERNAL F,  the  name  of   the   function   which   evaluates *
 *   the function to be minimized.                                           *
 *                                                                           *
 *   Input, INTEGER N, the number of variables.                              *
 *                                                                           *
 *   Input/output, DOUBLE PRECISION START(N). On  input,  a  starting  point *
 *   for the iteration. On output, this data may have been overwritten.      *
 *                                                                           *
 *   Output,  DOUBLE PRECISION XMIN(N),  the  coordinates  of   the    point *
 *   which is estimated to minimize the function.                            *
 *                                                                           *
 *   Output, DOUBLE PRECISION YNEWLO, the minimum value of the function.     *
 *                                                                           *
 *   Input, DOUBLE PRECISION REQMIN, the terminating limit for the  variance *
 *   of function values.                                                     *
 *                                                                           *
 *   Input,  DOUBLE PRECISION STEP(N),  determines  the  size   and    shape *
 *   of the  initial  simplex.  The  relative  magnitudes  of  its  elements *
 *   should reflect the units of the variables.                              *
 *                                                                           *
 *   Input,  INTEGER KONVGE,  the  convergence   check   is   carried    out *
 *   every KONVGE iterations.                                                *
 *                                                                           *
 *   Input, INTEGER KCOUNT, the maximum number of function evaluations.      *
 *                                                                           *
 *   Output, INTEGER ICOUNT, the number of function evaluations used.        *
 *                                                                           *
 *   Output, INTEGER NUMRES, the number of restarts.                         *
 *                                                                           *
 *   Output, INTEGER IFAULT, error indicator:                                *
 *     (0) No errors detected.                                               *
 *     (1) REQMIN, N, or KONVGE has an illegal value.                        *
 *     (2) Iteration  terminated  because  KCOUNT   was   exceeded   without *
 *         convergence.                                                      *
 *                                                                           *
 * Remark:                                                                   *
 *   The Nelder-Mead algorithm is also known as  "Downhill  simplex  method" *
 *   or "Amoeba method".                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* vim:set nu et ts=4 sw=4: */
