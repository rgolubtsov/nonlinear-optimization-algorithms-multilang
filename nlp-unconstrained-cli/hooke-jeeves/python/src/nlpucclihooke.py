# -*- coding: utf-8 -*-
# =============================================================================
# nlp-unconstrained-cli/hooke-jeeves/python/src/nlpucclihooke.py
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
# =============================================================================

from nlpucclihooko.rosenbrock import Rosenbrock
from nlpucclihooko.woods      import Woods

class NLPUCCLIHooke:
    """The NLPUCCLIHooke class contains methods to start up computations
    for solving a nonlinear optimization problem using the algorithm
    of Hooke and Jeeves with respect to two test problems:
    Rosenbrock and Woods.
    """

    def rosenbrock(self):
        """Solves the Rosenbrock test problem."""

        # Instantiating the Rosenbrock class.
        r = Rosenbrock()

        # Firing up computations.
        r.main()

        return None

    def woods(self):
        """Solves the Woods test problem."""

        # Instantiating the Woods class.
        w = Woods()

        # Firing up computations.
        w.main()

        return None

    def __init__(self):
        """Default constructor."""

        self = []

        return None

# =============================================================================
# vim:set nu:et:ts=4:sw=4:
# =============================================================================
