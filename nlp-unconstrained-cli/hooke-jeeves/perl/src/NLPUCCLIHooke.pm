#
# nlp-unconstrained-cli/hooke-jeeves/perl/src/NLPUCCLIHooke.pm
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
# =============================================================================
# Written by Radislav (Radicchio) Golubtsov, 2015-2023
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

=head1 NAME

NLPUCCLIHooke - The Hooke-Jeeves nonlinear minimization algorithm -
main application class

=head1 VERSION

Version 0.1.1

=head1 DESCRIPTION

The C<NLPUCCLIHooke> class contains methods to start up computations
for solving a nonlinear optimization problem using the algorithm
of Hooke and Jeeves with respect to two test problems:
Rosenbrock (see L<NLPUCCLIHooke::Rosenbrock>)
and Woods (see L<NLPUCCLIHooke::Woods>).

=cut

package NLPUCCLIHooke;

use strict;
use warnings;
use v5.10;

use NLPUCCLIHooke::Rosenbrock;
use NLPUCCLIHooke::Woods;

## The module version number.
our $VERSION = "0.1";

=head1 METHODS

=head2 rosenbrock

Solve the Rosenbrock test problem.

 @param {String[]} args - The array of command-line arguments.

=cut

sub rosenbrock {
    my $self = shift();
    my ($args) = @_;

    # Instantiating the Rosenbrock class.
    my $r = NLPUCCLIHooke::Rosenbrock->new();

    # Firing up computations.
    $r->main($args);

    return 1;
}

=head2 woods

Solve the Woods test problem.

 @param {String[]} args - The array of command-line arguments.

=cut

sub woods {
    my $self = shift();
    my ($args) = @_;

    # Instantiating the Woods class.
    my $w = NLPUCCLIHooke::Woods->new();

    # Firing up computations.
    $w->main($args);

    return 1;
}

1;

=head1 AUTHOR

Radislav (Radicchio) Golubtsov,
L<https://github.com/rgolubtsov|https://github.com/rgolubtsov>

=head1 COPYRIGHT AND LICENSE

There is no copyright - this code is in the public domain.

See the LICENSE or UNLICENSE files for more details in the root directory
of the project.

=cut

# vim:set nu et ts=4 sw=4:
