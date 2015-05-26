# =============================================================================
# nlp-unconstrained-core/hooke-jeeves/perl/src/NLPUCCoreHooke.pm
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
# =============================================================================

=head1 NAME

NLPUCCoreHooke - The Hooke-Jeeves nonlinear minimization algorithm -
main application class

=head1 VERSION

Version 0.1

=head1 DESCRIPTION

The C<NLPUCCoreHooke> class contains methods to start up computations
for solving a nonlinear optimization problem using the algorithm
of Hooke and Jeeves with respect to two test problems:
Rosenbrock (see L<NLPUCCoreHooke::Rosenbrock>)
and Woods (see L<NLPUCCoreHooke::Woods>).

=cut

package NLPUCCoreHooke;

use strict;
use warnings;
use v5.10;

use NLPUCCoreHooke::Rosenbrock;
use NLPUCCoreHooke::Woods;

## The app version number.
our $VERSION = "0.1";

=head1 ATTRIBUTES

=head2 FUN_EVALS

Global var. The number of function evaluations.

=cut

our $FUN_EVALS = 0;

=head1 METHODS

=head2 run_rosenbrock

Solve the Rosenbrock test problem.

 @param {String[]} args - The array of command-line arguments.

=cut

sub run_rosenbrock {
    my $self = shift();
    my ($args) = @_;

    # Instantiating the Rosenbrock class.
    my $rosenbrock = NLPUCCoreHooke::Rosenbrock->new();

    # Firing up computations.
    $rosenbrock->main($args);

    return 1;
}

=head2 run_woods

Solve the Woods test problem.

 @param {String[]} args - The array of command-line arguments.

=cut

sub run_woods {
    my $self = shift();
    my ($args) = @_;

    # Instantiating the Woods class.
    my $woods = NLPUCCoreHooke::Woods->new();

    # Firing up computations.
    $woods->main($args);

    return 1;
}

1;

=head1 AUTHOR

Radislav (Radic) Golubtsov,
L<https://github.com/rgolubtsov|https://github.com/rgolubtsov>

=head1 COPYRIGHT AND LICENSE

There is no copyright - this code is in the public domain.

See the LICENSE or UNLICENSE files for more details in the root directory
of the project.

=cut

# =============================================================================
# vim:set nu:et:ts=4:sw=4:
# =============================================================================
