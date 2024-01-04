#
# nlp-unconstrained-cli/hooke-jeeves/perl/src/NLPUCCLIHooke/Hooke.pm
# =============================================================================
# Nonlinear Optimization Algorithms Multilang. Version 0.1.1
# =============================================================================
# Nonlinear programming algorithms as the (un-)constrained minimization
# problems with the focus on their numerical expression using various
# programming languages.
#
# This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
# =============================================================================
# Written by Radislav (Radicchio) Golubtsov, 2015-2024
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

NLPUCCLIHooke::Hooke - The Hooke-Jeeves nonlinear minimization algorithm

=head1 DESCRIPTION

The C<Hooke> class contains methods for solving a nonlinear optimization
problem using the algorithm of Hooke and Jeeves.

=cut

package NLPUCCLIHooke::Hooke;

use strict;
use warnings;
use v5.10;

use Exporter "import";

## Helper constants.
use constant INDEX_ZERO      => 0;
use constant INDEX_ONE       => 1;
use constant ZERO_POINT_FIVE => 0.5;

=head1 ATTRIBUTES

=head2 VARS

Constant. The maximum number of variables.

=cut

use constant VARS => 250;

=head2 EPSMIN

Constant. The ending value of stepsize.

=cut

use constant EPSMIN => 1E-6;

=head2 IMAX

Constant. The maximum number of iterations.

=cut

use constant IMAX => 5000;

=head2 FUNEVALS

The number of function evaluations.

=cut

our $FUNEVALS = 0;

## Props to export.
our @EXPORT = (
    "INDEX_ZERO",
    "INDEX_ONE",
    "VARS",
    "EPSMIN",
    "IMAX",
);

=head1 METHODS

=head2 best_nearby

Helper method.

Given a point, look for a better one nearby, one coord at a time.

 @param   {Number[]} delta    - The delta between prevbest and point.
 @param   {Number[]} point    - The coordinate from where to begin.
 @param   {Number}   prevbest - The previous best-valued coordinate.
 @param   {Number}   nvars    - The number of variables.
 @param   {Object}   f        - The user-supplied objective function f(x,n).

 @returns {Number}              The objective function value at a nearby.

=cut

sub best_nearby {
    my ($delta, $point, $prevbest, $nvars, $f) = @_;

    my $minf;
    my @z = (0)x(VARS);
    my $ftmp;

    my $i;

    $minf = $prevbest;

    for ($i = 0; $i < $nvars; $i++) {
        $z[$i] = $point->[$i];
    }

    for ($i = 0; $i < $nvars; $i++) {
        $z[$i] = $point->[$i] + $delta->[$i];

        $ftmp = $f->(\@z, $nvars);

        if ($ftmp < $minf) {
            $minf = $ftmp;
        } else {
            $delta->[$i] = 0.0 - $delta->[$i];
            $z[$i]       = $point->[$i] + $delta->[$i];

            $ftmp = $f->(\@z, $nvars);

            if ($ftmp < $minf) {
                $minf = $ftmp;
            } else {
                $z[$i] = $point->[$i];
            }
        }
    }

    for ($i = 0; $i < $nvars; $i++) {
        $point->[$i] = $z[$i];
    }

    return $minf;
}

=head2 hooke

Main optimization method.

The hooke subroutine itself.

 @param   {Number}   nvars   - The number of variables.
 @param   {Number[]} startpt - The starting point coordinates.
 @param   {Number[]} endpt   - The ending point coordinates.
 @param   {Number}   rho     - The rho value.
 @param   {Number}   epsilon - The epsilon value.
 @param   {Number}   itermax - The maximum number of iterations.
 @param   {Object}   f       - The user-supplied objective function f(x,n).

 @returns {Number}             The number of iterations used
                               to find the local minimum.

=cut

sub hooke {
    my $self = shift();
    my ($nvars, $startpt, $endpt, $rho, $epsilon, $itermax, $f) = @_;

    my $i;
    my $iadj;
    my $iters;
    my $j;
    my $keep;

    my @newx    = (0)x(VARS);
    my @xbefore = (0)x(VARS);
    my @delta   = (0)x(VARS);
    my $steplength;
    my $fbefore;
    my $newf;
    my $tmp;

    for ($i = 0; $i < $nvars; $i++) {
        $xbefore[$i] = $startpt->[$i];
        $newx[$i]    = $xbefore[$i];

        $delta[$i] = abs($startpt->[$i] * $rho);

        if ($delta[$i] == 0.0) {
            $delta[$i] = $rho;
        }
    }

    $iadj       = 0;
    $steplength = $rho;
    $iters      = 0;

    $fbefore = $f->(\@newx, $nvars);

    $newf = $fbefore;

    while (($iters < $itermax) && ($steplength > $epsilon)) {
        $iters++;
        $iadj++;

        printf(
            "\nAfter %5d funevals, f(x) =  %.4le at\n", $FUNEVALS, $fbefore
        );

        for ($j = 0; $j < $nvars; $j++) {
            printf("   x[%2d] = %.4le\n", $j, $xbefore[$j]);
        }

        # Find best new point, one coord at a time.
        for ($i = 0; $i < $nvars; $i++) {
            $newx[$i] = $xbefore[$i];
        }

        $newf = best_nearby(\@delta, \@newx, $fbefore, $nvars, $f);

        # If we made some improvements, pursue that direction.
        $keep = 1;

        while (($newf < $fbefore) && ($keep == 1)) {
            $iadj = 0;

            for ($i = 0; $i < $nvars; $i++) {
                # Firstly, arrange the sign of delta[].
                if ($newx[$i] <= $xbefore[$i]) {
                    $delta[$i] = 0.0 - abs($delta[$i]);
                } else {
                    $delta[$i] = abs($delta[$i]);
                }

                # Now, move further in this direction.
                $tmp         = $xbefore[$i];
                $xbefore[$i] = $newx[$i];
                $newx[$i]    = $newx[$i] + $newx[$i] - $tmp;
            }

            $fbefore = $newf;

            $newf = best_nearby(\@delta, \@newx, $fbefore, $nvars, $f);

            # If the further (optimistic) move was bad....
            if ($newf >= $fbefore) {
                last;
            }

            # Make sure that the differences between the new and the old
            # points are due to actual displacements; beware of roundoff
            # errors that might cause newf < fbefore.
            $keep = 0;

            for ($i = 0; $i < $nvars; $i++) {
                $keep = 1;

                if (abs($newx[$i] - $xbefore[$i])
                    > (ZERO_POINT_FIVE * abs($delta[$i]))) {

                    last;
                } else {
                    $keep = 0;
                }
            }
        }

        if (($steplength >= $epsilon) && ($newf >= $fbefore)) {
            $steplength = $steplength * $rho;

            for ($i = 0; $i < $nvars; $i++) {
                $delta[$i] *= $rho;
            }
        }
    }

    for ($i = 0; $i < $nvars; $i++) {
        $endpt->[$i] = $xbefore[$i];
    }

    return $iters;
}

=head2 new

The Hooke class constructor.

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
