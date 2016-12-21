# -*- coding: utf-8 -*-
# =============================================================================
# nlp-unconstrained-core/hooke-jeeves/python/src/nlpuccorehooko/hooke.py
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
# =============================================================================

from nlpuccorehooko.funevals import FunEvals

class Hooke:
    """The Hooke class contains methods for solving a nonlinear optimization
    problem using the algorithm of Hooke and Jeeves.
    """

    ## Helper constants.
    INDEX_ZERO      = 0
    INDEX_ONE       = 1
    ZERO_POINT_FIVE = 0.5

    ## Constant. The maximum number of variables.
    VARS = 250

    ## Constant. The ending value of stepsize.
    EPSMIN = 1E-6

    ## Constant. The maximum number of iterations.
    IMAX = 5000

    def best_nearby(self, delta, point, prevbest, nvars, f, c_funevals):
        """Helper method.

        Given a point, look for a better one nearby, one coord at a time.

        Args:
            delta:      The delta between prevbest and point.
            point:      The coordinate from where to begin.
            prevbest:   The previous best-valued coordinate.
            nvars:      The number of variables.
            f:          The user-supplied objective function f(x,n).
            c_funevals: The number of function evaluations container
                        (FunEvals).

        Returns:
            The objective function value at a nearby.
        """

        z = [0]*self.VARS

        minf = prevbest

        i = 0

        while (i < nvars):
            z[i] = point[i]

            i += 1

        i = 0

        while (i < nvars):
            z[i] = point[i] + delta[i]

            ftmp = f(z, nvars, c_funevals)

            if (ftmp < minf):
                minf = ftmp
            else:
                delta[i] = 0.0 - delta[i]
                z[i]     = point[i] + delta[i]

                ftmp = f(z, nvars, c_funevals)

                if (ftmp < minf):
                    minf = ftmp
                else:
                    z[i] = point[i]

            i += 1

        i = 0

        while (i < nvars):
            point[i] = z[i]

            i += 1

        return minf

    def hooke(self, nvars, startpt, endpt, rho, epsilon, itermax, f):
        """Main optimization method.

        The hooke subroutine itself.

        Args:
            nvars:   The number of variables.
            startpt: The starting point coordinates.
            endpt:   The ending point coordinates.
            rho:     The rho value.
            epsilon: The epsilon value.
            itermax: The maximum number of iterations.
            f:       The user-supplied objective function f(x,n).

        Returns:
            The number of iterations used to find the local minimum.
        """

        newx    = [0]*self.VARS
        xbefore = [0]*self.VARS
        delta   = [0]*self.VARS

        i = 0

        while (i < nvars):
            xbefore[i] = startpt[i]
            newx[i]    = xbefore[i]

            delta[i] = abs(startpt[i] * rho)

            if (delta[i] == 0.0):
                delta[i] = rho

            i += 1

        iadj       = 0
        steplength = rho
        iters      = 0

        # Instantiating the FunEvals class.
        fe = FunEvals()

        fbefore = f(newx, nvars, fe)

        newf = fbefore

        while ((iters < itermax) and (steplength > epsilon)):
            iters += 1
            iadj  += 1

            print("\nAfter {:5d} funevals, f(x) =  {:.4e} at"
                .format(fe.get_funevals(), fbefore))

            j = 0

            while (j < nvars):
                print("   x[{:2d}] = {:.4e}".format(j, xbefore[j]))

                j += 1

            # Find best new point, one coord at a time.
            i = 0

            while (i < nvars):
                newx[i] = xbefore[i]

                i += 1

            newf = self.best_nearby(delta, newx, fbefore, nvars, f, fe)

            # If we made some improvements, pursue that direction.
            keep = 1

            while ((newf < fbefore) and (keep == 1)):
                iadj = 0

                i = 0

                while (i < nvars):
                    # Firstly, arrange the sign of delta[].
                    if (newx[i] <= xbefore[i]):
                        delta[i] = 0.0 - abs(delta[i])
                    else:
                        delta[i] = abs(delta[i])

                    # Now, move further in this direction.
                    tmp        = xbefore[i]
                    xbefore[i] = newx[i]
                    newx[i]    = newx[i] + newx[i] - tmp

                    i += 1

                fbefore = newf

                newf = self.best_nearby(delta, newx, fbefore, nvars, f, fe)

                # If the further (optimistic) move was bad....
                if (newf >= fbefore):
                    break

                # Make sure that the differences between the new and the old
                # points are due to actual displacements; beware of roundoff
                # errors that might cause newf < fbefore.
                keep = 0

                i = 0

                while (i < nvars):
                    keep = 1

                    if (abs(newx[i] - xbefore[i])
                        > (self.ZERO_POINT_FIVE * abs(delta[i]))):

                        break
                    else:
                        keep = 0

                    i += 1

            if ((steplength >= epsilon) and (newf >= fbefore)):
                steplength = steplength * rho

                i = 0

                while (i < nvars):
                    delta[i] *= rho

                    i += 1

        i = 0

        while (i < nvars):
            endpt[i] = xbefore[i]

            i += 1

        return iters

    def __init__(self):
        """Default constructor."""

        self = []

        return None

# =============================================================================
# vim:set nu:et:ts=4:sw=4:
# =============================================================================
