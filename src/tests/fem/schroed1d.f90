module schroed1d_assembly
use types, only: dp
use sorting, only: sort
use constants, only: pi
implicit none
private
public assemble_1d, Z, exact_energies, f

real(dp), parameter :: Z = 1._dp

contains

    pure function gaussian_potential(x, alpha, x0) result(V)
    ! Returns a Gaussian potential. See gaussian_density() for details.
    real(dp), intent(in) :: x(:), alpha, x0
    real(dp) :: V(size(x)), r(size(x))
    r = abs(x-x0)
    V = 2*alpha**2*r*erf(alpha*r) + 2*alpha*exp(-alpha**2*r**2)/sqrt(pi)
    end function

real(dp) elemental function f(x)
real(dp), intent(in) :: x
real(dp) :: tmp(1)
tmp = gaussian_potential([x], 2._dp, 4._dp)
f = tmp(1)
end function

subroutine assemble_1d(xin, nodes, ib, xiq, wtq, phihq, dphihq, Vq, Am, Bm)
! Assemble on a 2D rectangular uniform mesh
real(dp), intent(in):: xin(:), nodes(:), xiq(:), wtq(:), &
    phihq(:, :), dphihq(:, :), Vq(:,:)
integer, intent(in):: ib(:, :)
real(dp), intent(out):: Am(:,:), Bm(:, :)
real(dp), dimension(size(xiq), &
    size(xin)) :: phi_v, phi_dx, phi_dy
real(dp), dimension(size(xiq)) :: fq
real(dp), dimension(size(xiq)) :: x, y, xp, yp
integer :: Ne, p, e, i, j, iqx, iqy
real(dp) :: lx, ly
integer :: ax, ay, bx, by
real(dp) :: jacx, jacy, jac_det

Ne = size(nodes)-1
p = size(xin) - 1
! 2D shape functions
do ax = 1, p+1
    do iqx = 1, size(xiq)
        phi_v (iqx, ax) =  phihq(iqx, ax)
        phi_dx(iqx, ax) = dphihq(iqx, ax)
    end do
end do
Am=0; Bm=0
! Precalculate as much as possible:
lx = nodes(2) - nodes(1) ! Element size
jacx = lx/2
jac_det = abs(jacx)
xp = (xiq + 1) * jacx
phi_dx = phi_dx / jacx
do e = 1, Ne
    x = xp + nodes(e)
    do iqx = 1, size(xiq)
        fq(iqx) = f(x(iqx))
    end do
    do bx = 1, p+1
        j = ib(bx, e)
        if (j==0) cycle
        do ax = 1, p+1
            i = ib(ax, e)
            if (i == 0) cycle
            if (j > i) cycle
            Am(i,j) = Am(i,j) + sum(phi_dx(:, ax)*phi_dx(:, bx) &
                * jac_det * wtq) / 2
            Am(i,j) = Am(i,j) + sum(Vq(:,e) * &
                phi_v(:, ax)*phi_v(:, bx) &
                * jac_det * wtq)
            Bm(i,j) = Bm(i,j) + sum(( &
                phi_v(:, ax) * phi_v(:, bx) &
                * jac_det * wtq))
        end do
    end do
end do
do j = 1, size(Am, 2)
    do i = 1, j-1
        Am(i, j) = Am(j, i)
        Bm(i, j) = Bm(j, i)
    end do
end do
end subroutine

integer pure function calc_size(nmax) result(s)
integer, intent(in) :: nmax
integer :: n, m
s = 0
do n = 1, nmax
    do m = -n+1, n-1
        s = s + 1
    end do
end do
end function

subroutine exact_energies(Z, E)
real(dp), intent(in) :: Z
real(dp), intent(out) :: E(:)
real(dp) :: tmp(calc_size(size(E)))
integer :: n, m, idx
idx = 0
do n = 1, size(E)
    do m = -n+1, n-1
        idx = idx + 1
        tmp(idx) = -Z**2 / (2*(n-1._dp/2)**2)
    end do
end do
E = tmp(:size(E))
end subroutine

end module


!------------------------------------------------------------------------------

program schroed1d

use types, only: dp
use fe_mesh, only: cartesian_mesh_2d, define_connect_tensor_2d
use schroed1d_assembly, only: assemble_1d, Z, exact_energies, f
use feutils, only: get_parent_nodes, get_parent_quad_pts_wts, phih, dphih
use linalg, only: eigh
use mesh, only: meshexp
use feutils, only: define_connect
use splines, only: iixmin, spline3pars, poly3
implicit none

integer :: Nn, Ne
! xe(i) is the 'x' coordinate of the i-th mesh node
real(dp), allocatable :: xe(:)
integer :: Nq, p, Nb
real(dp), allocatable :: xin(:), xiq(:), wtq(:), A(:, :), B(:, :), c(:, :), &
    lam(:), phihq(:, :), dphihq(:, :), E_exact(:), x(:), Vq(:,:)
integer, allocatable :: ib(:, :), in(:, :)
real(dp) :: L, jacx
integer :: i, j, e, iqx

Ne = 20
p = 63
Nq = p+1
L = 8  ! The size of the box in atomic units

Nn = Ne + 1
allocate(xe(Nn))
xe = meshexp(0._dp, L, 1._dp, Ne) ! uniform mesh on [0, L]

print *, "Number of nodes:", Nn
print *, "Number of elements:", Ne

allocate(xin(Nq), x(Nq), Vq(Nq,Ne))
call get_parent_nodes(2, p, xin)
allocate(xiq(Nq), wtq(Nq))
call get_parent_quad_pts_wts(2, Nq, xiq, wtq)
allocate(phihq(size(xiq), size(xin)))
allocate(dphihq(size(xiq), size(xin)))
! Tabulate parent basis at quadrature points
forall(i=1:size(xiq), j=1:size(xin))  phihq(i, j) =  phih(xin, j, xiq(i))
forall(i=1:size(xiq), j=1:size(xin)) dphihq(i, j) = dphih(xin, j, xiq(i))

allocate(in(p+1,Ne),ib(p+1,Ne))
call define_connect(3,3,Ne,p,in,ib)

Nb = maxval(ib)
print *, "p =", p
print *, "DOFs =", Nb
allocate(A(Nb, Nb), B(Nb, Nb), c(Nb, Nb), lam(Nb))

call load_potential(xe, xiq, Vq)

print *, "Assembling..."
call assemble_1d(xin, xe, ib, xiq, wtq, phihq, dphihq, Vq, A, B)
print *, "Solving..."
call eigh(A, B, lam, c)
print *, "Eigenvalues:"
allocate(E_exact(20))
call exact_energies(Z, E_exact)
do i = 1, min(Nb, 20)
    print "(i4, f20.12)", i, lam(i)
end do

contains

subroutine load_potential(xe, xiq, Vq)
real(dp), intent(in) :: xe(:), xiq(:)
real(dp), intent(out) :: Vq(:,:)
real(dp), allocatable :: Xn(:), Vn(:), c(:,:)
real(dp) :: jacx, x(size(xiq))
integer :: u, n, e, ip
! Load the numerical potential
n = 1024
allocate(Xn(n), Vn(n))
open(newunit=u, file="../fft/sch1d_grid.txt", status="old")
read(u, *) Xn
read(u, *) Vn ! skip: atomic potential
read(u, *) Vn ! skip: density
read(u, *) Vn
close(u)
! Interpolate using cubic splines
allocate(c(0:4, n-1))
call spline3pars(Xn, Vn, [2, 2], [0._dp, 0._dp], c)
ip = 0
do e = 1, size(xe)-1
    jacx=(xe(e+1)-xe(e))/2;
    x = xe(e) + (xiq + 1) * jacx
    do iqx = 1, size(xiq)
        ip = iixmin(x(iqx), Xn, ip)
        Vq(iqx, e) = poly3(x(iqx), c(:, ip))
    end do
end do
end subroutine

end program