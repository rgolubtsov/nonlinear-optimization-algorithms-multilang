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
