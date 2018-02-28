/*
 * nlp-unconstrained-cli/hooke-jeeves/c/src/funevals.h
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

/**
 * The <code>funevals.h</code> header file is for helper purposes only.
 * It declares the number of objective function evaluations variable,
 * wrapped by the C structure (and corresponding accessor methods).
 *
 * @author  Radislav (Radicchio) Golubtsov
 * @version 0.1
 * @since   hooke-jeeves 0.1
 */

#ifndef __C__FUNEVALS_H
#define __C__FUNEVALS_H

/**
 * The structure to hold the number of function evaluations variable.
 * It is just to simulate a C++ semantic model where a class
 * encapsulates properties.
 */
struct fun_evals {
    /** The number of function evaluations. */
    unsigned int funevals;
};

/**
 * Getter for <code>fun_evals.funevals</code>.
 *
 * @return The number of function evaluations.
 */
extern unsigned int get_funevals(const struct fun_evals *);

/**
 * Setter for <code>fun_evals.funevals</code>.
 *
 * @param __fun_evals The number of function evaluations container (struct *).
 * @param __funevals  The number of function evaluations.
 */
extern void set_funevals(struct fun_evals *, const unsigned int);

#endif /* __C__FUNEVALS_H */

/* vim:set nu et ts=4 sw=4: */
