! =============================================================================
! nlp-unconstrained-core/hooke-jeeves/f95/rosenbrock.f95
! =============================================================================
! Nonlinear Optimization Algorithms Multilang. Version 0.1
! =============================================================================
! Nonlinear programming algorithms as the (un-)constrained minimization
! problems with the focus on their numerical expression using various
! programming languages.
!
! This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
! =============================================================================
! Copyright (C) 2015 Radislav (Radic) Golubtsov

! =============================================================================
! Nonlinear Optimization using the algorithm of Hooke and Jeeves
! 12 February 1994 author: Mark G. Johnson
! 
! Find a point X where the nonlinear function f(X) has a local  minimum.  X
! is an n-vector and f(X) is a scalar. In mathematical notation f:  R^n  ->
! R^1. The objective function f() is not required  to  be  continuous.  Nor
! does f() need to be differentiable. The program does not use  or  require
! derivatives of f().
! 
! The software user supplies three things: a subroutine that computes f(X),
! an initial "starting guess" of the minimum point X, and  values  for  the
! algorithm convergence parameters. Then the program searches for  a  local
! minimum, beginning from the  starting  guess,  using  the  Direct  Search
! algorithm of Hooke and Jeeves.
! 
! This C program is adapted from the Algol pseudocode found  in  "Algorithm
! 178: Direct Search" by Arthur F. Kaupe Jr., Communications  of  the  ACM,
! Vol 6. p.313 (June 1963). It includes the improvements suggested by  Bell
! and Pike (CACM v.9, p.684, Sept 1966) and  those  of  Tomlin  and  Smith,
! "Remark on Algorithm 178" (CACM v.12). The original paper, which I  don't
! recommend as highly as the one by A.  Kaupe,  is:  R.  Hooke  and  T.  A.
! Jeeves, "Direct Search Solution of Numerical and  Statistical  Problems",
! Journal of the ACM, Vol.8, April 1961, pp.212-229.
! 
! Calling sequence:
!   int hooke(nvars, startpt, endpt, rho, epsilon, itermax)
!     nvars    {an integer}           This is the number of  dimensions  in
!                                     the domain of f(). It is  the  number
!                                     of coordinates of the starting  point
!                                     (and the minimum point.)
!     startpt  {an array of doubles}  This is the  user-supplied  guess  at
!                                     the minimum.
!     endpt    {an array of doubles}  This is the  location  of  the  local
!                                     minimum, calculated by the program
!     rho      {a double}             This is a  user-supplied  convergence
!                                     parameter (more detail below),  which
!                                     should be set to a value between  0.0
!                                     and 1.0. Larger values  of  rho  give
!                                     greater probability of convergence on
!                                     highly nonlinear functions, at a cost
!                                     of more function evaluations. Smaller
!                                     values of rho reduces the  number  of
!                                     evaluations (and the program  running
!                                     time),  but  increases  the  risk  of
!                                     nonconvergence. See below.
!     epsilon  {a double}             This is the criterion for halting the
!                                     search  for  a  minimum.   When   the
!                                     algorithm begins  to  make  less  and
!                                     less progress on each  iteration,  it
!                                     checks the halting criterion: if  the
!                                     stepsize is below epsilon,  terminate
!                                     the iteration and return the  current
!                                     best estimate of the minimum.  Larger
!                                     values of epsilon  (such  as  1.0e-4)
!                                     give quicker running time, but a less
!                                     accurate  estimate  of  the  minimum.
!                                     Smaller values of  epsilon  (such  as
!                                     1.0e-7) give longer running time, but
!                                     a  more  accurate  estimate  of   the
!                                     minimum.
!     itermax  {an integer}           A  second,   rarely   used,   halting
!                                     criterion. If the algorithm  uses  >=
!                                     itermax iterations, halt.
! 
! The user-supplied objective function f(x,n) should return a  C  "double".
! Its arguments are x -- an array of doubles, and n -- an integer. x is the
! point at which  f(x)  should  be  evaluated,  and  n  is  the  number  of
! coordinates of x. That is, n is the number of coefficients being fitted.
! 
!   rho, the algorithm convergence control
! 
! The algorithm works by taking "steps" from one estimate of a minimum,  to
! another (hopefully better) estimate. Taking big steps gets to the minimum
! more quickly, at the risk of "stepping right over"  an  excellent  point.
! The stepsize is controlled by a user supplied parameter  called  rho.  At
! each iteration, the stepsize is multiplied by rho (0 < rho < 1),  so  the
! stepsize is successively reduced.
!     Small values of rho correspond to big stepsize  changes,  which  make
! the algorithm run more quickly. However, there is  a  chance  (especially
! with highly nonlinear functions) that these big changes will accidentally
! overlook a promising search vector, leading to nonconvergence.
!     Large values of rho correspond to small stepsize changes, which force
! the  algorithm  to   carefully   examine   nearby   points   instead   of
! optimistically  forging  ahead.  This   improves   the   probability   of
! convergence.
!     The stepsize is reduced until  it  is  equal  to  (or  smaller  than)
! epsilon. So  the  number  of  iterations  performed  by  Hooke-Jeeves  is
! determined by rho and epsilon:
! 
!   rho ** (number_of_iterations) = epsilon
! 
! In general it is a good idea to set rho to an  aggressively  small  value
! like 0.5 (hoping for fast convergence). Then, if the user  suspects  that
! the reported minimum is incorrect (or perhaps not accurate  enough),  the
! program can be run again with a larger value of rho such as  0.85,  using
! the result of the first minimization as the starting guess to  begin  the
! second minimization.
! 
! Normal use:
!   (1) Code your function f() in the C language
!   (2) Install your starting guess {or read it in}
!   (3) Run the program
!   (4) {for the skeptical}: Use the computed minimum as the starting point
!       for another run
! 
! Data Fitting:
!   Code your function f() to be the sum  of  the  squares  of  the  errors
!   (differences) between the computed values and the measured values. Then
!   minimize f() using Hooke-Jeeves.
!       EXAMPLE: you have 20 datapoints (ti, yi) and you want to find A,B,C
!   such that (A * t * t) + (B * exp(t)) + (C * tan(t)) fits  the  data  as
!   closely as possible. Then f() is just
! 
!   f(x) = SUM(measured_y[i] - ((A * t[i] * t[i]) +
!                               (B * exp(t[i])) +
!                               (C * tan(t[i])))) ^ 2
! 
!   where x[] is a 3-vector consisting of {A, B, C}.
! 
! The author of this software is M.G. Johnson.
! Permission to use, copy, modify, and distribute  this  software  for  any
! purpose without fee is hereby granted, provided that this  entire  notice
! is included in all copies of any software which is or includes a copy  or
! modification of this  software  and  in  all  copies  of  the  supporting
! documentation for such software. THIS SOFTWARE IS BEING PROVIDED "AS IS",
! WITHOUT ANY EXPRESS OR  IMPLIED  WARRANTY.  IN  PARTICULAR,  NEITHER  THE
! AUTHOR  NOR  AT&T  MAKE  ANY  REPRESENTATION  OR  WARRANTY  OF  ANY  KIND
! CONCERNING THE MERCHANTABILITY OF THIS SOFTWARE OR ITS  FITNESS  FOR  ANY
! PARTICULAR PURPOSE.
! =============================================================================

! Helper module.
module rosenbrock
    implicit none

    ! Max # of variables.
    integer, parameter :: VARS = 250

    ! Stepsize geometric shrink.
    double precision, parameter :: RHO_BEGIN = 0.5

    ! Ending value of stepsize.
    double precision, parameter :: EPSMIN = 1e-6

    ! Max # of iterations.
    integer, parameter :: IMAX = 5000

    ! Global variables.
    integer :: funevals = 0

    interface
        double precision function f(x, n)
            double precision, intent(in) :: x(:)

            integer, intent(inout), optional :: n
        end function f

        double precision function best_nearby(delta, point, prevbest, nvars)
            double precision, intent(inout) :: delta(:)
            double precision, intent(inout) :: point(:)
            double precision, intent(in)    :: prevbest

            integer, intent(in) :: nvars
        end function best_nearby

        integer function hooke(nvars, startpt, endpt, rho, epsilon, itermax)
            integer, intent(in) :: nvars

            double precision, intent(in)  :: startpt(:)
            double precision, intent(out) :: endpt(:)
            double precision, intent(in)  :: rho
            double precision, intent(in)  :: epsilon

            integer, intent(in) :: itermax
        end function hooke
    end interface
end module rosenbrock

! Rosenbrock's classic parabolic valley ("banana") function.
double precision function f(x, n)
    use rosenbrock, only: funevals

    implicit none

    double precision, intent(in) :: x(:)

    integer, intent(inout), optional :: n

    double precision :: a
    double precision :: b
    double precision :: c

    funevals = funevals + 1

    a = x(1)
    b = x(2)

    c = 100.0 * (b - (a * a)) * (b - (a * a))

    f = c + ((1.0 - a) * (1.0 - a))
end function f

! Given a point, look for a better one nearby, one coord at a time.
double precision function best_nearby(delta, point, prevbest, nvars)
    use rosenbrock, only: VARS, f

    implicit none

    double precision, intent(inout) :: delta(:)
    double precision, intent(inout) :: point(:)
    double precision, intent(in)    :: prevbest

    integer, intent(in) :: nvars

    double precision :: minf
    double precision :: z(VARS)
    double precision :: ftmp

    integer :: i

    minf = prevbest

    do i = 1, nvars
        z(i) = point(i)
    end do

    do i = 1, nvars
        z(i) = point(i) + delta(i)

        ftmp = f(z)

        if (ftmp < minf) then
            minf = ftmp
        else
            delta(i) = 0.0 - delta(i)
            z(i)     = point(i) + delta(i)

            ftmp = f(z)

            if (ftmp < minf) then
                minf = ftmp
            else
                z(i) = point(i)
            end if
        end if
    end do

    do i = 1, nvars
        point(i) = z(i)
    end do

    best_nearby = minf
end function best_nearby

integer function hooke(nvars, startpt, endpt, rho, epsilon, itermax)
    use rosenbrock, only: VARS, funevals, f, best_nearby

    implicit none

    integer, intent(in) :: nvars

    double precision, intent(in)  :: startpt(:)
    double precision, intent(out) :: endpt(:)
    double precision, intent(in)  :: rho
    double precision, intent(in)  :: epsilon

    integer, intent(in) :: itermax

    integer :: i
    integer :: iadj
    integer :: iters
    integer :: j
    integer :: keep

    double precision :: newx(VARS)
    double precision :: xbefore(VARS)
    double precision :: delta(VARS)
    double precision :: steplength
    double precision :: fbefore
    double precision :: newf
    double precision :: tmp

    do i = 1, nvars
        xbefore(i) = startpt(i)
        newx(i)    = xbefore(i)

        delta(i) = abs(startpt(i) * rho)

        if (delta(i) == 0.0) then
            delta(i) = rho
        end if
    end do

    iadj       = 0
    steplength = rho
    iters      = 0

    fbefore = f(newx)

    newf = fbefore

    do while ((iters < itermax) .and. (steplength > epsilon))
        iters = iters + 1
        iadj  = iadj + 1

        write ( &
            *, '(/, "After ", i5, " funevals, f(x) =  ", 1pe11.4e3, " at")' &
        ) funevals, fbefore

        do j = 1, nvars
            write (*, '("   x[", i2, "] = ", 1pe12.4e3)') j - 1, xbefore(j)
        end do

        ! Find best new point, one coord at a time.
        do i = 1, nvars
            newx(i) = xbefore(i)
        end do

        newf = best_nearby(delta, newx, fbefore, nvars)

        ! If we made some improvements, pursue that direction.
        keep = 1

        do while ((newf < fbefore) .and. (keep == 1))
            iadj = 0

            do i = 1, nvars
                ! Firstly, arrange the sign of delta[].
                if (newx(i) <= xbefore(i)) then
                    delta(i) = 0.0 - abs(delta(i))
                else
                    delta(i) = abs(delta(i))
                end if

                ! Now, move further in this direction.
                tmp        = xbefore(i)
                xbefore(i) = newx(i)
                newx(i)    = newx(i) + newx(i) - tmp
            end do

            fbefore = newf

            newf = best_nearby(delta, newx, fbefore, nvars)

            ! If the further (optimistic) move was bad....
            if (newf >= fbefore) then
                exit
            end if

            ! Make sure that the differences between the new and the old
            ! points are due to actual displacements; beware of roundoff
            ! errors that might cause newf < fbefore.
            keep = 0

            do i = 1, nvars
                keep = 1

                if (abs(newx(i) - xbefore(i)) > (0.5 * abs(delta(i)))) then
                    exit
                else
                    keep = 0
                end if
            end do
        end do

        if ((steplength >= epsilon) .and. (newf >= fbefore)) then
            steplength = steplength * rho

            do i = 1, nvars
                delta(i) = delta(i) * rho
            end do
        end if
    end do

    do i = 1, nvars
        endpt(i) = xbefore(i)
    end do

    hooke = iters
end function hooke

program banana
    use rosenbrock, only: VARS, RHO_BEGIN, EPSMIN, IMAX, hooke

    implicit none

    integer :: nvars
    integer :: itermax
    integer :: jj
    integer :: i

    double precision :: startpt(VARS)
    double precision :: rho
    double precision :: epsilon
    double precision :: endpt(VARS)

    ! Starting guess for Rosenbrock's test function.
    nvars      = 2
    startpt(1) = -1.2
    startpt(2) = 1.0
    itermax    = IMAX
    rho        = RHO_BEGIN
    epsilon    = EPSMIN

    jj = hooke(nvars, startpt, endpt, rho, epsilon, itermax)

    write (*, '(///, "HOOKE USED ", i2, " ITERATIONS, AND RETURNED")') jj

    do i = 1, nvars
        write (*, '("x[", i3, "] = ", 1pe15.7e3, " ")') i - 1, endpt(i)
    end do
end program banana

! =============================================================================
! vim:set nu:et:ts=4:sw=4:
! =============================================================================
