!
! nlp-unconstrained-cli/hooke-jeeves/f95/src/hooke.f95
! =============================================================================
! Nonlinear Optimization Algorithms Multilang. Version 0.1.1
! =============================================================================
! Nonlinear programming algorithms as the (un-)constrained minimization
! problems with the focus on their numerical expression using various
! programming languages.
!
! This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
! =============================================================================
! Written by Radislav (Radicchio) Golubtsov, 2015-2024
!
! This is free and unencumbered software released into the public domain.
!
! Anyone is free to copy, modify, publish, use, compile, sell, or
! distribute this software, either in source code form or as a compiled
! binary, for any purpose, commercial or non-commercial, and by any
! means.
!
! (See the LICENSE file at the top of the source tree.)
!

! === Helper module.
! Contains constants and variables declaration and function prototypes.
module proto__
    implicit none

    ! Constant. The double precision kind for real type values.
    integer, parameter :: DP = kind(1.d0)

    ! Constant. The maximum number of variables.
    integer, parameter :: VARS = 250

    ! Constant. The stepsize geometric shrink.
    real(DP), parameter :: RHO_BEGIN = 0.5

    ! The Hooke & Jeeves algorithm works reasonably well on Rosenbrock's
    ! function, but can fare worse on some standard test functions,
    ! depending on rho. Here is an example that works well when rho = 0.5,
    ! but fares poorly with rho = 0.6, and better again with rho = 0.8.
    real(DP), parameter :: RHO_WOODS = 0.6

    ! Constant. The ending value of stepsize.
    real(DP), parameter :: EPSMIN = 1e-6

    ! Constant. The maximum number of iterations.
    integer, parameter :: IMAX = 5000

    ! The number of function evaluations.
    integer :: FUNEVALS = 0

    ! Function prototypes.
    interface
        ! The objective function f(x,n).
        real(kind(1.d0)) function f(x, n)
            real(kind(1.d0)), intent(in) :: x(:)
            integer, intent(in)          :: n
        end function f

        ! The helper function best_nearby(...).
        real(kind(1.d0)) function best_nearby(delta, point, prevbest, nvars)
            real(kind(1.d0)), intent(inout) :: delta(:)
            real(kind(1.d0)), intent(inout) :: point(:)
            real(kind(1.d0)), intent(in)    :: prevbest
            integer, intent(in)             :: nvars
        end function best_nearby

        ! The main optimization function hooke(...).
        integer function hooke(nvars, startpt, endpt, rho, epsilon, itermax)
            integer, intent(in)           :: nvars
            real(kind(1.d0)), intent(in)  :: startpt(:)
            real(kind(1.d0)), intent(out) :: endpt(:)
            real(kind(1.d0)), intent(in)  :: rho
            real(kind(1.d0)), intent(in)  :: epsilon
            integer, intent(in)           :: itermax
        end function hooke
    end interface
end module proto__

! === The user-supplied objective function f(x,n).
real(DP) function f(x, n)
    use proto__, only: DP, FUNEVALS

    implicit none

    ! Returns: The objective function value.

    ! Arg. The point at which f(x) should be evaluated.
    real(DP), intent(in) :: x(:)

    ! Arg. The number of coordinates of x.
    integer, intent(in) :: n

#ifndef WOODS
    ! Rosenbrock's classic parabolic valley ("banana") function.
    real(DP) :: a
    real(DP) :: b
    real(DP) :: c

    FUNEVALS = FUNEVALS + 1

    a = x(1)
    b = x(2)

    c = 100.0 * (b - (a * a)) * (b - (a * a))

    f = c + ((1.0 - a) * (1.0 - a))
#else
    ! Woods -- a la More, Garbow & Hillstrom (TOMS algorithm 566).
    real(DP) :: s1
    real(DP) :: s2
    real(DP) :: s3
    real(DP) :: t1
    real(DP) :: t2
    real(DP) :: t3
    real(DP) :: t4
    real(DP) :: t5

    FUNEVALS = FUNEVALS + 1

    s1 = x(2) - x(1) * x(1)
    s2 = 1    - x(1)
    s3 = x(2) - 1

    t1 = x(4) - x(3) * x(3)
    t2 = 1    - x(3)
    t3 = x(4) - 1

    t4 = s3 + t3
    t5 = s3 - t3

    f = 100 * (s1 * s1) + s2 * s2 &
       + 90 * (t1 * t1) + t2 * t2 &
       + 10 * (t4 * t4) + t5 * t5 / 10.
#endif

    if (n == 0) then
        write (*, '("Warning: The number of coordinates of x = 0.")')
    end if
end function f

! === Helper function.
! Given a point, look for a better one nearby, one coord at a time.
real(DP) function best_nearby(delta, point, prevbest, nvars)
    use proto__, only: DP, VARS, f

    implicit none

    ! Returns: The objective function value at a nearby.

    ! Arg. The delta between prevbest and point.
    real(DP), intent(inout) :: delta(:)

    ! Arg. The coordinate from where to begin.
    real(DP), intent(inout) :: point(:)

    ! Arg. The previous best-valued coordinate.
    real(DP), intent(in) :: prevbest

    ! Arg. The number of variables.
    integer, intent(in) :: nvars

    real(DP) :: minf
    real(DP) :: z(VARS)
    real(DP) :: ftmp

    integer :: i

    minf = prevbest

    do i = 1, nvars
        z(i) = point(i)
    end do

    do i = 1, nvars
        z(i) = point(i) + delta(i)

        ftmp = f(z, nvars)

        if (ftmp < minf) then
            minf = ftmp
        else
            delta(i) = 0.0 - delta(i)
            z(i)     = point(i) + delta(i)

            ftmp = f(z, nvars)

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

! === Main optimization function.
! The hooke subroutine itself.
integer function hooke(nvars, startpt, endpt, rho, epsilon, itermax)
    use proto__, only: DP, VARS, FUNEVALS, f, best_nearby

    implicit none

    ! Returns: The number of iterations used to find the local minimum.

    ! Arg. The number of variables.
    integer, intent(in) :: nvars

    ! Arg. The starting point coordinates.
    real(DP), intent(in) :: startpt(:)

    ! Arg. The ending point coordinates.
    real(DP), intent(out) :: endpt(:)

    ! Arg. The rho value.
    real(DP), intent(in) :: rho

    ! Arg. The epsilon value.
    real(DP), intent(in) :: epsilon

    ! Arg. The maximum number of iterations.
    integer, intent(in) :: itermax

    integer :: i
    integer :: iadj
    integer :: iters
    integer :: j
    integer :: keep

    real(DP) :: newx(VARS)
    real(DP) :: xbefore(VARS)
    real(DP) :: delta(VARS)
    real(DP) :: steplength
    real(DP) :: fbefore
    real(DP) :: newf
    real(DP) :: tmp

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

    fbefore = f(newx, nvars)

    newf = fbefore

    do while ((iters < itermax) .and. (steplength > epsilon))
        iters = iters + 1
        iadj  = iadj + 1

        write ( &
            *, '(/, "After ", i5, " funevals, f(x) =  ", 1pe10.4e2, " at")' &
        ) FUNEVALS, fbefore

        do j = 1, nvars
            write (*, '("   x[", i2, "] = ", 1pe11.4e2)') j - 1, xbefore(j)
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
                ! Firstly, arrange the sign of delta().
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

! === Main program.
program hooke__
    use proto__, only: DP, VARS, RHO_BEGIN, RHO_WOODS, EPSMIN, IMAX, hooke

    implicit none

    integer :: nvars
    integer :: itermax
    integer :: jj
    integer :: i

    real(DP) :: startpt(VARS)
    real(DP) :: rho
    real(DP) :: epsilon
    real(DP) :: endpt(VARS)

#ifndef WOODS
    ! Starting guess for Rosenbrock's test function.
    nvars      =  2
    startpt(1) = -1.2
    startpt(2) =  1.0
    rho        =  RHO_BEGIN
#else
    ! Starting guess test problem "Woods".
    nvars      =  4
    startpt(1) = -3
    startpt(2) = -1
    startpt(3) = -3
    startpt(4) = -1
    rho        =  RHO_WOODS
#endif

    itermax = IMAX
    epsilon = EPSMIN

    jj = hooke(nvars, startpt, endpt, rho, epsilon, itermax)

    write (*, '(///, "HOOKE USED ", i2, " ITERATIONS, AND RETURNED")') jj

    do i = 1, nvars
        write (*, '("x[", i3, "] = ", 1pe15.7e2, " ")') i - 1, endpt(i)
    end do

#ifdef WOODS
    write (*, '("True answer: f(1, 1, 1, 1) = 0.")')
#endif
end program hooke__

! vim:set nu et ts=4 sw=4:
