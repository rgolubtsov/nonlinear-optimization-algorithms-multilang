# -*- coding: utf-8 -*-
# =============================================================================
# nlp-unconstrained-core/hooke-jeeves/python/src/nlpuccorehooko/woods.py
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
# =============================================================================

from nlpuccorehooko.hooke import Hooke

class Woods:
    """The Woods class is responsible for solving a nonlinear optimization
    problem using the algorithm of Hooke and Jeeves.

    The objective function in this case is the so-called Woods function.
    """

    ## Helper constants.
    INDEX_TWO   =   2
    INDEX_THREE =   3
    ONE_HUNDRED = 100
    NINETY      =  90
    TEN         =  10
    TEN_POINT   =  10.
    FOUR        =   4
    MINUS_THREE =  -3
    MINUS_ONE   =  -1

    ##
    # Constant. The stepsize geometric shrink.
    #
    # The Hooke & Jeeves algorithm works reasonably well on Rosenbrock's
    # function, but can fare worse on some standard test functions,
    # depending on rho. Here is an example that works well when rho = 0.5,
    # but fares poorly with rho = 0.6, and better again with rho = 0.8.
    #
    RHO_WOODS = 0.6

    def f(self, x, n):
        """The user-supplied objective function f(x,n).

        Woods -- a la More, Garbow & Hillstrom (TOMS algorithm 566).

        Args:
            x: The point at which f(x) should be evaluated.
            n: The number of coordinates of x.

        Returns:
            The objective function value.
        """

        Hooke.funEvals += 1

        s1 = x[Hooke.INDEX_ONE]  - x[Hooke.INDEX_ZERO] * x[Hooke.INDEX_ZERO]
        s2 = 1                   - x[Hooke.INDEX_ZERO]
        s3 = x[Hooke.INDEX_ONE]  - 1

        t1 = x[self.INDEX_THREE] - x[self.INDEX_TWO]   * x[self.INDEX_TWO]
        t2 = 1                   - x[self.INDEX_TWO]
        t3 = x[self.INDEX_THREE] - 1

        t4 = s3 + t3
        t5 = s3 - t3

        return (self.ONE_HUNDRED * (s1 * s1) + s2 * s2
                   + self.NINETY * (t1 * t1) + t2 * t2
                      + self.TEN * (t4 * t4) + t5 * t5 / self.TEN_POINT)

    def main(self):
        """Main program function main() :-)."""

        startpt = [0]*Hooke.VARS
        endpt   = [0]*Hooke.VARS

        # Starting guess test problem "Woods".
        nvars                      =  self.FOUR
        startpt[Hooke.INDEX_ZERO]  =  self.MINUS_THREE
        startpt[Hooke.INDEX_ONE]   =  self.MINUS_ONE
        startpt[ self.INDEX_TWO]   =  self.MINUS_THREE
        startpt[ self.INDEX_THREE] =  self.MINUS_ONE
        itermax                    = Hooke.IMAX
        rho                        =  self.RHO_WOODS
        epsilon                    = Hooke.EPSMIN

        # Instantiating the Hooke class.
        h = Hooke()

        jj = h.hooke(nvars, startpt, endpt, rho, epsilon, itermax, self.f)

        print("\n\n\nHOOKE USED " + str(jj) + " ITERATIONS, AND RETURNED")

        i = 0

        while (i < nvars):
            print("x[{:3d}] = {:15.7e} ".format(i, endpt[i]))

            i += 1

        print("True answer: f(1, 1, 1, 1) = 0.")

        return None

    def __init__(self):
        """Default constructor."""

        self = []

        return None

# =============================================================================
# vim:set nu:et:ts=4:sw=4:
# =============================================================================
