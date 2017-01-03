/*
 * nlp-unconstrained-cli/hooke-jeeves/vala/src/funevals.vala
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

namespace CLIHooke {

/**
 * The <code>FunEvals</code> class is a helper class.
 * It holds the only property &ndash; the number of objective function
 * evaluations (and corresponding accessor methods).
 */
class FunEvals {
    /** The number of function evaluations. */
    uint funevals;

    /**
     * Getter for <code>funevals</code>.
     *
     * @return The number of function evaluations.
     */
    public uint get_funevals() {
        return funevals;
    }

    /**
     * Setter for <code>funevals</code>.
     *
     * @param __funevals The number of function evaluations.
     */
    public void set_funevals(uint __funevals) {
        funevals = __funevals;
    }

    /** Default constructor. */
    public FunEvals() {}
}

} // namespace CLIHooke

// vim:set nu:et:ts=4:sw=4:
