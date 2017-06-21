module mixings

! This module contains SCF mixing algorithms.

use types, only: dp
use utils, only: stop_error
implicit none
private
public mixing_linear, mixing_anderson

interface
    subroutine F_fn(x, y, energies)
    ! y = F(x), also return the calculated energies to converge
    import :: dp
    implicit none
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:), energies(:)
    end subroutine

    real(dp) function integral_fn(x)
    ! Computes the integral of the vector 'x'
    import :: dp
    implicit none
    real(dp), intent(in) :: x(:)
    end function
end interface

contains

subroutine mixing_linear(F, integral, x0, nenergies, max_iter, alpha, eps, &
        x_out)
! Finds "x" so that F(x) = x
procedure(F_fn) :: F
procedure(integral_fn) :: integral
real(dp), intent(in) :: x0(:)
integer, intent(in) :: nenergies, max_iter
real(dp), intent(in) :: alpha
real(dp), intent(in) :: eps
real(dp), intent(out) :: x_out(:)

real(dp), dimension(size(x0)) :: x_i, y_i, R_i
real(dp) :: old_energies(nenergies), energies(nenergies)
real(dp) :: x_i_norm, R_i_norm
real(dp) :: err_old, err, L2_err
integer :: i
x_i = x0
err_old = 1e12_dp
old_energies = 1e12_dp
do i = 1, max_iter
    call F(x_i, y_i, energies)
    R_i = y_i-x_i

    ! L2 norm of the "input" potential:
    x_i_norm = sqrt(integral(x_i**2))
    ! L2 norm of the "output-input" potential:
    R_i_norm = sqrt(integral(R_i**2))
    if (x_i_norm < 1e-12_dp) x_i_norm = 1e-12_dp
    L2_err = R_i_norm / x_i_norm
    err = maxval(abs(energies - old_energies))
    ! Do at least 3 iterations
    if (i >= 3 .and. L2_err < 5e-5_dp) then
        if (err < eps .and. err_old < eps) then
            x_out = x_i
            return
        end if
    end if
    old_energies = energies
    err_old = err

    x_i = x_i + alpha * R_i
end do
call stop_error("SCF didn't converge")
end subroutine

subroutine mixing_anderson(F, integral, x0, nenergies, max_iter, alpha, eps, &
        x_out)
! Finds "x" so that F(x) = x, uses x0 as the initial estimate
procedure(F_fn) :: F
procedure(integral_fn) :: integral
real(dp), intent(in) :: x0(:)
integer, intent(in) :: nenergies, max_iter
real(dp), intent(in) :: alpha
real(dp), intent(in) :: eps
real(dp), intent(out) :: x_out(:)

real(dp), dimension(size(x0)) :: x_i, y_i, x1_i, R_i, R1_i, delta_R, delta_x
real(dp) :: beta
real(dp) :: sn, sd
real(dp) :: old_energies(nenergies), energies(nenergies)
real(dp) :: x_i_norm, R_i_norm
real(dp) :: err_old, err, L2_err
integer :: i
x_i = x0
err_old = 1e12_dp
old_energies = 1e12_dp
do i = 1, max_iter
    call F(x_i, y_i, energies)
    R_i = y_i-x_i

    ! L2 norm of the "input" potential:
    x_i_norm = sqrt(integral(x_i**2))
    ! L2 norm of the "output-input" potential:
    R_i_norm = sqrt(integral(R_i**2))
    if (x_i_norm < 1e-12_dp) x_i_norm = 1e-12_dp
    L2_err = R_i_norm / x_i_norm
    err = maxval(abs(energies - old_energies))
    ! Do at least 3 iterations
    if (i >= 3 .and. L2_err < 5e-5_dp) then
        if (err < eps .and. err_old < eps) then
            x_out = x_i
            return
        end if
    end if
    old_energies = energies
    err_old = err

    if (i > 1) then
        delta_x = x_i - x1_i
        delta_R = R_i - R1_i
    end if
    x1_i = x_i
    R1_i = R_i
    x_i = x_i + alpha * R_i
    if (i > 1) then
        sn = integral(R_i * delta_R)
        sd = integral(delta_R**2)
        beta = sn / sd
        x_i = x_i - beta * (delta_x + alpha * delta_R)
    end if
end do
call stop_error("SCF didn't converge")
end subroutine


end module
