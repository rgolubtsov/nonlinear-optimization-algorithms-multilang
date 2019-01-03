#
# nlp-unconstrained-cli/hooke-jeeves/perl/src/NLPUCCLIHooke/Rosenbrock.pm
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
# =============================================================================
# Written by Radislav (Radicchio) Golubtsov, 2015-2019
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

NLPUCCLIHooke::Rosenbrock - The Hooke-Jeeves nonlinear minimization algorithm
in application to the Rosenbrock test problem

=head1 DESCRIPTION

The C<Rosenbrock> class is responsible for solving a nonlinear optimization
problem using the algorithm of Hooke and Jeeves (Rosenbrock test problem).

The objective function in this case is the Rosenbrock's parabolic valley
function.

=cut

package NLPUCCLIHooke::Rosenbrock;

use strict;
use warnings;
use v5.10;

use NLPUCCLIHooke::Hooke;

## Helper constants.
use constant ONE_HUNDRED_POINT_ZERO =>  100.0;
use constant ONE_POINT_ZERO         =>  1.0;
use constant TWO                    =>  2;
use constant MINUS_ONE_POINT_TWO    => -1.2;

=head1 ATTRIBUTES

=head2 RHO_BEGIN

Constant. The stepsize geometric shrink.

=cut

use constant RHO_BEGIN => 0.5;

=head1 METHODS

=head2 f

The user-supplied objective function f(x,n).

Represents here the Rosenbrock's classic parabolic valley ("banana") function.

 @param   {Number[]} x - The point at which f(x) should be evaluated.
 @param   {Number}   n - The number of coordinates of x.

 @returns {Number}       The objective function value.

=cut

sub f {
    my ($x, $n) = @_;

    my $a;
    my $b;
    my $c;

    $NLPUCCLIHooke::Hooke::FUNEVALS++;

    $a = $x->[INDEX_ZERO];
    $b = $x->[INDEX_ONE];

    $c = ONE_HUNDRED_POINT_ZERO * ($b - ($a * $a)) * ($b - ($a * $a));

    return ($c + ((ONE_POINT_ZERO - $a) * (ONE_POINT_ZERO - $a)));
}

=head2 main

Main program function.

 @param {String[]} args - The array of command-line arguments.

=cut

sub main {
    my $self = shift();
    my ($args) = @_;

    my $nvars;
    my $itermax;
    my $jj;
    my $i;

    my @startpt = (0)x(VARS);
    my $rho;
    my $epsilon;
    my @endpt   = (0)x(VARS);

    # Starting guess for Rosenbrock's test function.
    $nvars               = TWO;
    $startpt[INDEX_ZERO] = MINUS_ONE_POINT_TWO;
    $startpt[INDEX_ONE]  = ONE_POINT_ZERO;
    $itermax             = IMAX;
    $rho                 = RHO_BEGIN;
    $epsilon             = EPSMIN;

    # Instantiating the Hooke class.
    my $h = NLPUCCLIHooke::Hooke->new();

    $jj = $h->hooke($nvars, \@startpt, \@endpt, $rho, $epsilon, $itermax, \&f);

    say("\n\n\nHOOKE USED $jj ITERATIONS, AND RETURNED");

    for ($i = 0; $i < $nvars; $i++) {
        printf("x[%3d] = %15.7le \n", $i, $endpt[$i]);
    }

    return 1;
}

=head2 new

The Rosenbrock class constructor.

 @returns {Object} The object instance of the class.

=cut

sub new {
    my $class = shift();
    my $self  = [];

    bless($self, $class);

    return $self;
}

1;

# vim:set nu et ts=4 sw=4:
