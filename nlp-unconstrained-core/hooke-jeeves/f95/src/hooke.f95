! =============================================================================
! nlp-unconstrained-core/hooke-jeeves/f95/src/hooke.f95
! =============================================================================
! Nonlinear Optimization Algorithms Multilang. Version 0.1
! =============================================================================
! Nonlinear programming algorithms as the (un-)constrained minimization
! problems with the focus on their numerical expression using various
! programming languages.
!
! This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
! =============================================================================

! === Helper module.
! Contains constants and variables declaration and function prototypes.
module proto__
    implicit none

    ! Constant. The maximum number of variables.
    integer, parameter :: VARS = 250

    ! Constant. The stepsize geometric shrink.
    double precision, parameter :: RHO_BEGIN = 0.5

    ! The Hooke & Jeeves algorithm works reasonably well on Rosenbrock's
    ! function, but can fare worse on some standard test functions,
    ! depending on rho. Here is an example that works well when rho = 0.5,
    ! but fares poorly with rho = 0.6, and better again with rho = 0.8.
    double precision, parameter :: RHO_WOODS = 0.6

    ! Constant. The ending value of stepsize.
    double precision, parameter :: EPSMIN = 1e-6

    ! Constant. The maximum number of iterations.
    integer, parameter :: IMAX = 5000

    ! The number of function evaluations.
    integer :: funevals = 0

    ! Function prototypes.
    interface
        ! The objective function f(x,n).
        double precision function f(x, n)
            double precision, intent(in)     :: x(:)
            integer, intent(inout), optional :: n
        end function f

        ! The helper function best_nearby(...).
        double precision function best_nearby(delta, point, prevbest, nvars)
            double precision, intent(inout) :: delta(:)
            double precision, intent(inout) :: point(:)
            double precision, intent(in)    :: prevbest
            integer, intent(in)             :: nvars
        end function best_nearby

        ! The main optimization function hooke(...).
        integer function hooke(nvars, startpt, endpt, rho, epsilon, itermax)
            integer, intent(in)           :: nvars
            double precision, intent(in)  :: startpt(:)
            double precision, intent(out) :: endpt(:)
            double precision, intent(in)  :: rho
            double precision, intent(in)  :: epsilon
            integer, intent(in)           :: itermax
        end function hooke
    end interface
end module proto__

! === The user-supplied objective function f(x,n).
double precision function f(x, n)
    use proto__, only: funevals

    implicit none

    ! Returns: The objective function value.

    ! Arg. The point at which f(x) should be evaluated.
    double precision, intent(in) :: x(:)

    ! Arg. The number of coordinates of x.
    integer, intent(inout), optional :: n

#ifndef WOODS
! Rosenbrock's classic parabolic valley ("banana") function.
    double precision :: a
    double precision :: b
    double precision :: c

    funevals = funevals + 1

    a = x(1)
    b = x(2)

    c = 100.0 * (b - (a * a)) * (b - (a * a))

    f = c + ((1.0 - a) * (1.0 - a))
#else
! Woods -- a la More, Garbow & Hillstrom (TOMS algorithm 566).
    double precision :: s1
    double precision :: s2
    double precision :: s3
    double precision :: t1
    double precision :: t2
    double precision :: t3
    double precision :: t4
    double precision :: t5

    funevals = funevals + 1

    s1 = x(2) - x(1) * x(1)
    s2 = 1 - x(1)
    s3 = x(2) - 1
    t1 = x(4) - x(3) * x(3)
    t2 = 1 - x(3)
    t3 = x(4) - 1
    t4 = s3 + t3
    t5 = s3 - t3

    f = 100 * (s1 * s1) + s2 * s2 &
       + 90 * (t1 * t1) + t2 * t2 &
       + 10 * (t4 * t4) + t5 * t5 / 10.
#endif
end function f

! === Helper function.
! Given a point, look for a better one nearby, one coord at a time.
double precision function best_nearby(delta, point, prevbest, nvars)
    use proto__, only: VARS, f

    implicit none

    ! Returns: The objective function value at a nearby.

    ! Arg. The delta between prevbest and point.
    double precision, intent(inout) :: delta(:)

    ! Arg. The coordinate from where to begin.
    double precision, intent(inout) :: point(:)

    ! Arg. The previous best-valued coordinate.
    double precision, intent(in) :: prevbest

    ! Arg. The number of variables.
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

! === Main optimization function.
! The hooke subroutine itself.
integer function hooke(nvars, startpt, endpt, rho, epsilon, itermax)
    use proto__, only: VARS, funevals, f, best_nearby

    implicit none

    ! Returns: The number of iterations used to find the local minimum.

    ! Arg. The number of variables.
    integer, intent(in) :: nvars

    ! Arg. The starting point coordinates.
    double precision, intent(in) :: startpt(:)

    ! Arg. The ending point coordinates.
    double precision, intent(out) :: endpt(:)

    ! Arg. The rho value.
    double precision, intent(in) :: rho

    ! Arg. The epsilon value.
    double precision, intent(in) :: epsilon

    ! Arg. The maximum number of iterations.
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
    use proto__, only: VARS, RHO_BEGIN, RHO_WOODS, EPSMIN, IMAX, hooke

    implicit none

    integer :: nvars
    integer :: itermax
    integer :: jj
    integer :: i

    double precision :: startpt(VARS)
    double precision :: rho
    double precision :: epsilon
    double precision :: endpt(VARS)

#ifndef WOODS
    ! Starting guess for Rosenbrock's test function.
    nvars      = 2
    startpt(1) = -1.2
    startpt(2) = 1.0
    rho        = RHO_BEGIN
#else
    ! Starting guess test problem "Woods".
    nvars      = 4
    startpt(1) = -3
    startpt(2) = -1
    startpt(3) = -3
    startpt(4) = -1
    rho        = RHO_WOODS
#endif

    itermax    = IMAX
    epsilon    = EPSMIN

    jj = hooke(nvars, startpt, endpt, rho, epsilon, itermax)

    write (*, '(///, "HOOKE USED ", i2, " ITERATIONS, AND RETURNED")') jj

    do i = 1, nvars
        write (*, '("x[", i3, "] = ", 1pe15.7e3, " ")') i - 1, endpt(i)
    end do

#ifdef WOODS
    write (*, '("True answer: f(1, 1, 1, 1) = 0.")')
#endif
end program hooke__

! =============================================================================
! vim:set nu:et:ts=4:sw=4:
! =============================================================================
