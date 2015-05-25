#!/usr/bin/env perl
# =============================================================================
# nlp-unconstrained-core/hooke-jeeves/perl/bin/hooke
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1
# =============================================================================

use strict;
use warnings;
use v5.10;

use NLPUCCoreHooke;

# Solving the Rosenbrock test problem.
NLPUCCoreHooke->run_rosenbrock();

# =============================================================================
# vim:set nu:et:ts=4:sw=4:
# =============================================================================
