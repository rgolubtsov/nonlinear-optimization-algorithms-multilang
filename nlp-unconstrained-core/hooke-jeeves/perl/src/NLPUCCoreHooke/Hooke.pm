# =============================================================================
# nlp-unconstrained-core/hooke-jeeves/perl/src/NLPUCCoreHooke/Hooke.pm
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

NLPUCCoreHooke::Hooke - The Hooke-Jeeves nonlinear minimization algorithm

=head1 DESCRIPTION

The C<Hooke> class contains methods for solving a nonlinear optimization
problem using the algorithm of Hooke and Jeeves.

=cut

package NLPUCCoreHooke::Hooke;

use strict;
use warnings;
use v5.10;

use Exporter "import";

## Helper constant.
use constant INDEX_ZERO => 0;

## Helper constant.
use constant INDEX_ONE => 1;

## Helper constant.
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

 @param   {Number[]} delta     - The delta between prev_best and point.
 @param   {Number[]} point     - The coordinate from where to begin.
 @param   {Number}   prev_best - The previous best-valued coordinate.
 @param   {Number}   n_vars    - The number of variables.
 @param   {Object}   f         - The user-supplied objective function f(x,n).

 @returns {Number}               The objective function value at a nearby.

=cut

sub best_nearby {
    my ($delta, $point, $prev_best, $n_vars, $f) = @_;

    my $min_f;
    my @z = (0)x(VARS);
    my $f_tmp;

    my $i;

    $min_f = $prev_best;

    for ($i = 0; $i < $n_vars; $i++) {
        $z[$i] = $point->[$i];
    }

    for ($i = 0; $i < $n_vars; $i++) {
        $z[$i] = $point->[$i] + $delta->[$i];

        $f_tmp = $f->(\@z, $n_vars);

        if ($f_tmp < $min_f) {
            $min_f = $f_tmp;
        } else {
            $delta->[$i] = 0.0 - $delta->[$i];
            $z[$i]       = $point->[$i] + $delta->[$i];

            $f_tmp = $f->(\@z, $n_vars);

            if ($f_tmp < $min_f) {
                $min_f = $f_tmp;
            } else {
                $z[$i] = $point->[$i];
            }
        }
    }

    for ($i = 0; $i < $n_vars; $i++) {
        $point->[$i] = $z[$i];
    }

    return $min_f;
}

=head2 hooke

Main optimization method.

The hooke subroutine itself.

 @param   {Number}   n_vars   - The number of variables.
 @param   {Number[]} start_pt - The starting point coordinates.
 @param   {Number[]} end_pt   - The ending point coordinates.
 @param   {Number}   rho      - The rho value.
 @param   {Number}   epsilon  - The epsilon value.
 @param   {Number}   iter_max - The maximum number of iterations.
 @param   {Object}   f        - The user-supplied objective function f(x,n).

 @returns {Number}              The number of iterations used
                                to find the local minimum.

=cut

sub hooke {
    my $self = shift();
    my ($n_vars, $start_pt, $end_pt, $rho, $epsilon, $iter_max, $f) = @_;

    my $i;
    my $i_adj;
    my $iters;
    my $j;
    my $keep;

    my @new_x    = (0)x(VARS);
    my @x_before = (0)x(VARS);
    my @delta    = (0)x(VARS);
    my $step_length;
    my $f_before;
    my $new_f;
    my $tmp;

    for ($i = 0; $i < $n_vars; $i++) {
        $x_before[$i] = $start_pt->[$i];
        $new_x[$i]    = $x_before[$i];

        $delta[$i] = abs($start_pt->[$i] * $rho);

        if ($delta[$i] == 0.0) {
            $delta[$i] = $rho;
        }
    }

    $i_adj       = 0;
    $step_length = $rho;
    $iters       = 0;

    $f_before = $f->(\@new_x, $n_vars);

    $new_f = $f_before;

    while (($iters < $iter_max) && ($step_length > $epsilon)) {
        $iters++;
        $i_adj++;

        printf(
            "\nAfter %5d funevals, f(x) =  %.4le at\n",
            $NLPUCCoreHooke::FUN_EVALS, $f_before
        );

        for ($j = 0; $j < $n_vars; $j++) {
            printf("   x[%2d] = %.4le\n", $j, $x_before[$j]);
        }

        # Find best new point, one coord at a time.
        for ($i = 0; $i < $n_vars; $i++) {
            $new_x[$i] = $x_before[$i];
        }

        $new_f = best_nearby(\@delta, \@new_x, $f_before, $n_vars, $f);

        # If we made some improvements, pursue that direction.
        $keep = 1;

        while (($new_f < $f_before) && ($keep == 1)) {
            $i_adj = 0;

            for ($i = 0; $i < $n_vars; $i++) {
                # Firstly, arrange the sign of delta[].
                if ($new_x[$i] <= $x_before[$i]) {
                    $delta[$i] = 0.0 - abs($delta[$i]);
                } else {
                    $delta[$i] = abs($delta[$i]);
                }

                # Now, move further in this direction.
                $tmp          = $x_before[$i];
                $x_before[$i] = $new_x[$i];
                $new_x[$i]    = $new_x[$i] + $new_x[$i] - $tmp;
            }

            $f_before = $new_f;

            $new_f = best_nearby(\@delta, \@new_x, $f_before, $n_vars, $f);

            # If the further (optimistic) move was bad....
            if ($new_f >= $f_before) {
                last;
            }

            # Make sure that the differences between the new and the old
            # points are due to actual displacements; beware of roundoff
            # errors that might cause new_f < f_before.
            $keep = 0;

            for ($i = 0; $i < $n_vars; $i++) {
                $keep = 1;

                if (abs($new_x[$i] - $x_before[$i])
                    > (ZERO_POINT_FIVE * abs($delta[$i]))) {

                    last;
                } else {
                    $keep = 0;
                }
            }
        }

        if (($step_length >= $epsilon) && ($new_f >= $f_before)) {
            $step_length = $step_length * $rho;

            for ($i = 0; $i < $n_vars; $i++) {
                $delta[$i] *= $rho;
            }
        }
    }

    for ($i = 0; $i < $n_vars; $i++) {
        $end_pt->[$i] = $x_before[$i];
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

# =============================================================================
# vim:set nu:et:ts=4:sw=4:
# =============================================================================
