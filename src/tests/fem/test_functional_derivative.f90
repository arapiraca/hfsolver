program test_functional_derivative

! Test for 4 Gaussian charges (based off of test_free_energy4.f90) and a
! functional derivative propagation.

use types, only: dp
use ofdft_fe, only: free_energy2, fe_data, initialize_fe, &
    free_energy2_low_level, free_energy
use ofdft_fft, only: reciprocal_space_vectors, radial_potential_fourier, &
    real2fourier
use constants, only: Ha2eV, pi, i_
use utils, only: loadtxt, assert, linspace, zeros
use splines, only: spline3pars, iixmin, poly3
use isolve, only: solve_cg
use interp3d, only: trilinear
use feutils, only: quad_lobatto
use md, only: positions_fcc
use converged_energies, only: four_gaussians
use poisson3d_assembly, only: func2quad, integral, assemble_3d_coo_rhs
use fe_mesh, only: c2fullc_3d, fe2quad_3d, fe2quad_3d_lobatto
implicit none
real(dp) :: Eee, Een, Ts, Exc, Etot, Etot_conv
integer :: p, DOF, Nq
real(dp) :: Rcut, L, T_eV, T_au
integer, parameter :: natom = 4
real(dp) :: X(3, natom)
integer :: Nex, Ney, Nez

real(dp), allocatable, dimension(:, :, :, :) :: nenq_pos, nq_pos
type(fe_data) :: fed

real(dp), allocatable, dimension(:, :, :, :) :: Hpsi, &
    psi, Venq, &
    nenq_neutral
complex(dp), allocatable, dimension(:, :, :, :) :: cpsi, cpsi2, cpsi3
real(dp) :: dt
integer :: max_iter
real(dp) :: brent_eps, free_energy_
real(dp) :: psi_norm
real(dp) :: background
real(dp), allocatable :: rhs(:), sol(:), fullsol(:), fullsol2(:)
real(dp) :: Nelec
integer :: i

Rcut = 0.3_dp
p = 8
Nex = 8
Ney = 8
Nez = 8
L = 2
T_eV = 0.0862_dp
T_au = T_ev / Ha2eV
Nq = 9
call positions_fcc(X, L)
call initialize_fe(L, Nex, Ney, Nez, p, Nq, quad_lobatto, fed)

allocate(nenq_pos(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(nq_pos(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
nenq_pos = func2quad(fed%nodes, fed%elems, fed%xiq, nen)
nq_pos = func2quad(fed%nodes, fed%elems, fed%xiq, fne)

Nelec = real(natom, dp)

brent_eps = 1e-3_dp
max_iter = 200

allocate(nenq_neutral(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(Venq(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(Hpsi(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(psi(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(cpsi(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(cpsi2(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(cpsi3(fed%Nq, fed%Nq, fed%Nq, fed%Ne))


! Calculate Venq
allocate(rhs(fed%Nb), sol(fed%Nb), fullsol(maxval(fed%in)))
allocate(fullsol2(maxval(fed%in)))
background = integral(fed%nodes, fed%elems, fed%wtq3, nenq_pos) / &
    (fed%Lx*fed%Ly*fed%Lz)
print *, "Total (negative) ionic charge: ", background * (fed%Lx*fed%Ly*fed%Lz)
print *, "Subtracting constant background (Q/V): ", background
nenq_neutral = nenq_pos - background
print *, "Assembling RHS..."
call assemble_3d_coo_rhs(fed%Ne, fed%p, 4*pi*nenq_neutral, fed%jac_det, fed%wtq3, &
    fed%ib, fed%phi_v, rhs)
print *, "sum(rhs):    ", sum(rhs)
print *, "integral rhs:", integral(fed%nodes, fed%elems, fed%wtq3, nenq_neutral)
print *, "Solving..."
sol = solve_cg(fed%Ap, fed%Aj, fed%Ax, rhs, zeros(size(rhs)), 1e-12_dp, 800)
print *, "Converting..."
call c2fullc_3d(fed%in, fed%ib, sol, fullsol)
if (fed%spectral) then
    call fe2quad_3d_lobatto(fed%elems, fed%xiq, fed%in, fullsol, Venq)
else
    call fe2quad_3d(fed%elems, fed%xin, fed%xiq, fed%phihq, fed%in, fullsol, Venq)
end if
print *, "Done"


psi = sqrt(nq_pos)
psi_norm = integral(fed%nodes, fed%elems, fed%wtq3, psi**2)
print *, "Initial norm of psi:", psi_norm
psi = sqrt(Nelec / psi_norm) * psi
psi_norm = integral(fed%nodes, fed%elems, fed%wtq3, psi**2)
print *, "norm of psi:", psi_norm
! This returns H[n] = delta F / delta n, we save it to the Hpsi variable to
! save space:
call free_energy(fed%nodes, fed%elems, fed%in, fed%ib, fed%Nb, fed%Lx, fed%Ly, fed%Lz, fed%xin, fed%xiq, fed%wtq3, T_au, &
    Venq, psi**2, fed%phihq, fed%Ap, fed%Aj, fed%Ax, fed%matd, fed%spectral, &
    fed%phi_v, fed%jac_det, &
    Eee, Een, Ts, Exc, free_energy_, Hpsi=Hpsi)
! Hpsi = H[psi] = delta F / delta psi = 2*H[n]*psi, due to d/dpsi = 2 psi d/dn
Hpsi = Hpsi * 2*psi

DOF = fed%Nb

Etot = Ts + Een + Eee + Exc
Etot_conv = sum(four_gaussians)
print *, "p =", p
print *, "DOF =", DOF
print *, "Rcut =", Rcut
print *, "T_au =", T_au
print *, "Summary of energies [a.u.]:"
print "('    Ts   = ', f14.8)", Ts
print "('    Een  = ', f14.8)", Een
print "('    Eee  = ', f14.8)", Eee
print "('    Exc  = ', f14.8)", Exc
print *, "   ---------------------"
print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV

print *, "Errors:"
print *, abs(Ts - four_gaussians(1))
print *, abs(Een - four_gaussians(2))
print *, abs(Eee - four_gaussians(3))
print *, abs(Exc - four_gaussians(4))
print *, abs(Etot - Etot_conv)
call assert(abs(Ts - four_gaussians(1)) < 1e-8_dp)
call assert(abs(Een - four_gaussians(2)) < 1e-8_dp)
call assert(abs(Eee - four_gaussians(3)) < 1e-8_dp)
call assert(abs(Exc - four_gaussians(4)) < 1e-8_dp)
call assert(abs(Etot - Etot_conv) < 1e-8_dp)

! Propagate

cpsi = psi
dt = 0.1_dp

! Do first step by hand:
print *, "First step"
cpsi2 = cpsi
cpsi = cpsi2 - i_*dt*Hpsi*cpsi2

do i = 1, 3
    print *, "iter =", i
    cpsi3 = cpsi2; cpsi2 = cpsi
    cpsi = cpsi3 - 2*i_*dt*Hpsi*cpsi2
    psi = abs(cpsi)
end do
print *, "Done"

contains

real(dp) function nen(x_, y_, z_) result(n)
real(dp), intent(in) :: x_, y_, z_
real(dp), parameter :: alpha = 6, Z = 1
real(dp) :: r2
integer :: i, a, b, c
n = 0
do i = 1, natom
    do a = -1, 1
    do b = -1, 1
    do c = -1, 1
        r2 = sum(([x_, y_, z_]-X(:, i)+[a, b, c]*L)**2)
        n = n - Z*alpha**3/pi**(3._dp/2)*exp(-alpha**2*r2)
    end do
    end do
    end do
end do
end function

real(dp) function fne(x, y, z) result(n)
real(dp), intent(in) :: x, y, z
real(dp), parameter :: alpha = 5, Z_ = 4
real(dp) :: r
r = sqrt((x-L/2)**2+(y-L/2)**2+(z-L/2)**2)
n = Z_*alpha**3/pi**(3._dp/2)*exp(-alpha**2*R**2)
end function

end program
