# -*- coding: utf-8 -*-
# nlp-unconstrained-cli/hooke-jeeves/python/src/nlpucclihooko/rosenbrock.py
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
# =============================================================================
# Written by Radislav (Radicchio) Golubtsov, 2015-2020
#
# This is free and unencumbered software released into the public domain.
#
# Anyone is free to copy, modify, publish, use, compile, sell, or
# distribute this software, either in source code form or as a compiled
# binary, for any purpose, commercial or non-commercial, and by any
# means.
#
# (See the LICENSE file at the top of the source tree.)
#

from nlpucclihooko.hooke import Hooke

class Rosenbrock:
    """The Rosenbrock class is responsible for solving a nonlinear optimization
    problem using the algorithm of Hooke and Jeeves.

    The objective function in this case is the Rosenbrock's parabolic valley
    function.
    """

    ## Helper constants.
    ONE_HUNDRED_POINT_ZERO = 100.0
    ONE_POINT_ZERO         =   1.0
    TWO                    =   2
    MINUS_ONE_POINT_TWO    =  -1.2

    ## Constant. The stepsize geometric shrink.
    RHO_BEGIN = 0.5

    def f(self, x, n, c_funevals):
        """The user-supplied objective function f(x,n).

        Represents here the Rosenbrock's classic parabolic valley
        ("banana") function.

        Args:
            x:          The point at which f(x) should be evaluated.
            n:          The number of coordinates of x.
            c_funevals: The number of function evaluations container
                        (FunEvals).

        Returns:
            The objective function value.
        """

        c_funevals.set_funevals(c_funevals.get_funevals() + 1)

        a = x[Hooke.INDEX_ZERO]
        b = x[Hooke.INDEX_ONE]

        c = self.ONE_HUNDRED_POINT_ZERO * (b - (a * a)) * (b - (a * a))

        return (c + ((self.ONE_POINT_ZERO - a) * (self.ONE_POINT_ZERO - a)))

    def main(self):
        """Main program function main() :-)."""

        startpt = [0]*Hooke.VARS
        endpt   = [0]*Hooke.VARS

        # Starting guess for Rosenbrock's test function.
        nvars                     =  self.TWO
        startpt[Hooke.INDEX_ZERO] =  self.MINUS_ONE_POINT_TWO
        startpt[Hooke.INDEX_ONE]  =  self.ONE_POINT_ZERO
        itermax                   = Hooke.IMAX
        rho                       =  self.RHO_BEGIN
        epsilon                   = Hooke.EPSMIN

        # Instantiating the Hooke class.
        h = Hooke()

        jj = h.hooke(nvars, startpt, endpt, rho, epsilon, itermax, self.f)

        print("\n\n\nHOOKE USED " + str(jj) + " ITERATIONS, AND RETURNED")

        i = 0

        while (i < nvars):
            print("x[{:3d}] = {:15.7e} ".format(i, endpt[i]))

            i += 1

        return None

    def __init__(self):
        """Default constructor."""

        self = []

        return None

# vim:set nu et ts=4 sw=4:
