! =============================================================================
! nlp-unconstrained-core/hooke-jeeves/f95/woods.f95
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

! Helper module.
module woods
    implicit none

    ! Max # of variables.
    integer, parameter :: VARS = 250

    ! The Hooke & Jeeves algorithm works reasonably well on Rosenbrock's
    ! function, but can fare worse on some standard test functions,
    ! depending on rho. Here is an example that works well when rho = 0.5,
    ! but fares poorly with rho = 0.6, and better again with rho = 0.8.
    double precision, parameter :: RHO_WOODS = 0.6

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
            double precision, intent(in) :: prevbest
            integer, intent(in) :: nvars
        end function best_nearby

        integer function hooke(nvars, startpt, endpt, rho, epsilon, itermax)
            integer, intent(in) :: nvars
            double precision, intent(in) :: startpt(:)
            double precision, intent(out) :: endpt(:)
            double precision, intent(in) :: rho
            double precision, intent(in) :: epsilon
            integer, intent(in) :: itermax
        end function hooke
    end interface
end module woods

! Woods -- a la More, Garbow & Hillstrom (TOMS algorithm 566).
double precision function f(x, n)
    use woods, only: funevals
    implicit none

    double precision, intent(in) :: x(:)
    integer, intent(inout), optional :: n

    double precision s1
    double precision s2
    double precision s3
    double precision t1
    double precision t2
    double precision t3
    double precision t4
    double precision t5

    funevals = funevals + 1
    s1 = x(2) - x(1) * x(1)
    s2 = 1 - x(1)
    s3 = x(2) - 1
    t1 = x(4) - x(3) * x(3)
    t2 = 1 - x(3)
    t3 = x(4) - 1
    t4 = s3 + t3
    t5 = s3 - t3

    f = 100 * (s1 * s1) + s2 * s2 + 90 * (t1 * t1) + t2 * t2 + 10 * (t4 * t4) &
        + t5 * t5 / 10.
end function f

! Given a point, look for a better one nearby, one coord at a time.
double precision function best_nearby(delta, point, prevbest, nvars)
    use woods, only: VARS, f
    implicit none

    double precision, intent(inout) :: delta(:)
    double precision, intent(inout) :: point(:)
    double precision, intent(in) :: prevbest
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
            z(i) = point(i) + delta(i)
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
    use woods, only: VARS, funevals, f, best_nearby
    implicit none

    integer, intent(in) :: nvars
    double precision, intent(in) :: startpt(:)
    double precision, intent(out) :: endpt(:)
    double precision, intent(in) :: rho
    double precision, intent(in) :: epsilon
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
        newx(i) = xbefore(i)
        delta(i) = abs(startpt(i) * rho)

        if (delta(i) == 0.0) then
            delta(i) = rho
        end if
    end do

    iadj = 0
    steplength = rho
    iters = 0
    fbefore = f(newx)
    newf = fbefore

    do while ((iters < itermax) .and. (steplength > epsilon))
        iters = iters + 1
        iadj = iadj + 1
        write (*, &
            '(/, "After ", i5, " funevals, f(x) =  ", 1pe11.4e3, " at")') &
            funevals, fbefore

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
                tmp = xbefore(i)
                xbefore(i) = newx(i)
                newx(i) = newx(i) + newx(i) - tmp
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

program more_garbow_hillstrom
    use woods, only: VARS, RHO_WOODS, EPSMIN, IMAX, hooke
    implicit none

    integer :: nvars
    integer :: itermax
    integer :: jj
    integer :: i
    double precision :: startpt(VARS)
    double precision :: rho
    double precision :: epsilon
    double precision :: endpt(VARS)

    ! Starting guess test problem "Woods".
    nvars = 4
    startpt(1) = -3
    startpt(2) = -1
    startpt(3) = -3
    startpt(4) = -1
    itermax = IMAX
    rho = RHO_WOODS
    epsilon = EPSMIN
    jj = hooke(nvars, startpt, endpt, rho, epsilon, itermax)
    write (*, '(///, "HOOKE USED ", i2, " ITERATIONS, AND RETURNED")') jj

    do i = 1, nvars
        write (*, '("x[", i3, "] = ", 1pe15.7e3, " ")') i - 1, endpt(i)
    end do

    write (*, '("True answer: f(1, 1, 1, 1) = 0.")')
end program more_garbow_hillstrom

! =============================================================================
! vim:set nu:et:ts=4:sw=4:
! =============================================================================
