# -*- coding: utf-8 -*-
# nlp-unconstrained-cli/hooke-jeeves/python/src/nlpucclihooke.py
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
# =============================================================================
# Written by Radislav (Radicchio) Golubtsov, 2015-2025
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

# vim:set nu et ts=4 sw=4:
