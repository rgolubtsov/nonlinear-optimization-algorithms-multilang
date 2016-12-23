#
# nlp-unconstrained-cli/hooke-jeeves/perl/src/NLPUCCLIHooke/Woods.pm
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
# =============================================================================
# Written by Radislav (Radicchio) Golubtsov, 2016
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

NLPUCCLIHooke::Woods - The Hooke-Jeeves nonlinear minimization algorithm
in application to the Woods test problem

=head1 DESCRIPTION

The C<Woods> class is responsible for solving a nonlinear optimization problem
using the algorithm of Hooke and Jeeves (Woods test problem).

The objective function in this case is the so-called I<Woods> function.

=cut

package NLPUCCLIHooke::Woods;

use strict;
use warnings;
use v5.10;

use NLPUCCLIHooke::Hooke;

## Helper constants.
use constant INDEX_TWO   =>  2;
use constant INDEX_THREE =>  3;
use constant ONE_HUNDRED =>  100;
use constant NINETY      =>  90;
use constant TEN         =>  10;
use constant TEN_POINT   =>  10.;
use constant FOUR        =>  4;
use constant MINUS_THREE => -3;
use constant MINUS_ONE   => -1;

=head1 ATTRIBUTES

=head2 RHO_WOODS

Constant. The stepsize geometric shrink.

 The Hooke & Jeeves algorithm works reasonably well on Rosenbrock's
 function, but can fare worse on some standard test functions,
 depending on rho. Here is an example that works well when rho = 0.5,
 but fares poorly with rho = 0.6, and better again with rho = 0.8.

=cut

use constant RHO_WOODS => 0.6;

=head1 METHODS

=head2 f

The user-supplied objective function f(x,n).

Woods -- a la More, Garbow & Hillstrom (TOMS algorithm 566).

 @param   {Number[]} x - The point at which f(x) should be evaluated.
 @param   {Number}   n - The number of coordinates of x.

 @returns {Number}       The objective function value.

=cut

sub f {
    my ($x, $n) = @_;

    my $s1;
    my $s2;
    my $s3;
    my $t1;
    my $t2;
    my $t3;
    my $t4;
    my $t5;

    $NLPUCCLIHooke::Hooke::FUNEVALS++;

    $s1 = $x->[INDEX_ONE]   - $x->[INDEX_ZERO] * $x->[INDEX_ZERO];
    $s2 = 1                 - $x->[INDEX_ZERO];
    $s3 = $x->[INDEX_ONE]   - 1;

    $t1 = $x->[INDEX_THREE] - $x->[INDEX_TWO]  * $x->[INDEX_TWO];
    $t2 = 1                 - $x->[INDEX_TWO];
    $t3 = $x->[INDEX_THREE] - 1;

    $t4 = $s3 + $t3;
    $t5 = $s3 - $t3;

    return (ONE_HUNDRED * ($s1 * $s1) + $s2 * $s2
               + NINETY * ($t1 * $t1) + $t2 * $t2
                  + TEN * ($t4 * $t4) + $t5 * $t5 / TEN_POINT);
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

    # Starting guess test problem "Woods".
    $nvars                = FOUR;
    $startpt[INDEX_ZERO]  = MINUS_THREE;
    $startpt[INDEX_ONE]   = MINUS_ONE;
    $startpt[INDEX_TWO]   = MINUS_THREE;
    $startpt[INDEX_THREE] = MINUS_ONE;
    $itermax              = IMAX;
    $rho                  = RHO_WOODS;
    $epsilon              = EPSMIN;

    # Instantiating the Hooke class.
    my $h = NLPUCCLIHooke::Hooke->new();

    $jj = $h->hooke($nvars, \@startpt, \@endpt, $rho, $epsilon, $itermax, \&f);

    say("\n\n\nHOOKE USED $jj ITERATIONS, AND RETURNED");

    for ($i = 0; $i < $nvars; $i++) {
        printf("x[%3d] = %15.7le \n", $i, $endpt[$i]);
    }

    say("True answer: f(1, 1, 1, 1) = 0.");

    return 1;
}

=head2 new

The Woods class constructor.

 @returns {Object} The object instance of the class.

=cut

sub new {
    my $class = shift();
    my $self  = [];

    bless($self, $class);

    return $self;
}

1;

# vim:set nu:et:ts=4:sw=4:
