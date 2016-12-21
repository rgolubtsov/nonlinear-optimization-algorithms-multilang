/*
 * nlp-unconstrained-cli/hooke-jeeves/cc/src/funevals.h
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

#ifndef __CC__FUNEVALS_H
#define __CC__FUNEVALS_H

/**
 * The <code>NLPUCCLIHooke</code> namespace is used as a container
 * for the <code>FunEvals</code> class.
 */
namespace NLPUCCLIHooke {

/**
 * The <code>FunEvals</code> class is a helper class.
 * It holds the only property &ndash; the number of objective function
 * evaluations (and corresponding accessor methods).
 *
 * @author  Radislav (Radic) Golubtsov
 * @version 0.1
 * @see     Hooke
 * @since   hooke-jeeves 0.1
 */
class FunEvals {
private:
    /** The number of function evaluations. */
    unsigned int funEvals;

public:
    /**
     * Getter for <code>funEvals</code>.
     *
     * @return The number of function evaluations.
     */
    unsigned int getFunEvals();

    /**
     * Setter for <code>funEvals</code>.
     *
     * @param __funEvals The number of function evaluations.
     */
    void setFunEvals(const unsigned int);

    /** Default constructor. */
    FunEvals();

    /** Destructor. */
    ~FunEvals();
};

} // namespace NLPUCCLIHooke

#endif // __CC__FUNEVALS_H

// ============================================================================
// vim:set nu:et:ts=4:sw=4:
// ============================================================================
