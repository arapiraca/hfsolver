module pksdft_fft
use constants, only: pi
use types, only: dp
use utils, only: assert, stop_error, clock, allocate_mold
use pofdft_fft, only: preal2fourier, pfourier2real
use arpack, only: peig, mul_matrix_vector
implicit none
private
public solve_schroedinger, applyH, solve_schroedinger2

contains

subroutine solve_schroedinger(myid, comm_all, commy, commz, Ng, nsub, Vloc, &
        L, G2, cutfn, nev, ncv, eigs, orbitals)
integer, intent(in) :: myid, comm_all, commy, commz, Ng(3), nsub(3), nev, ncv
real(dp), intent(in) :: Vloc(:,:,:) ! Local effective potential
real(dp), intent(in) :: G2(:,:,:), L(:), cutfn(:,:,:)
real(dp), intent(out) :: eigs(:) ! eigs(nev)
! orbitals(Ng_local(1),Ng_local(2),Ng_local(3),nev)
real(dp), intent(out) :: orbitals(:,:,:,:)

real(dp), allocatable :: d(:), v(:,:)
integer :: Ng_local(3), n
real(dp) :: t1, t2
logical :: verbose
verbose = .false.
Ng_local = Ng / nsub
n = product(Ng_local)
allocate(v(n,ncv), d(ncv))
call cpu_time(t1)
call peig(comm_all, myid, n, nev, ncv, "SA", av, d, v)
call cpu_time(t2)
eigs = d(:nev)+10
orbitals = reshape(v(:,:nev), [Ng_local(1),Ng_local(2),Ng_local(3),nev]) &
    * sqrt(product(Ng/L))
if (myid == 0 .and. verbose) then
    print *, "Arpack Time:", t2-t1
end if

contains

    subroutine av(x, y)
    ! Compute y = A*x
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    real(dp), dimension(Ng_local(1),Ng_local(2),Ng_local(3)) :: &
        psi
    psi = reshape(x, [Ng_local(1),Ng_local(2),Ng_local(3)])
    call applyH(commy, commz, Ng, nsub, Vloc-10, G2, cutfn, psi)
    y = reshape(psi, [product(Ng_local)])
    end

end subroutine


subroutine chebyshev_filter(X, Y, m, a, b, HX)
! Filter vectors in X by an `m` degree Chebyshev polynomial that dampens on the
! interval [a,b], and the output filtered vectors in Y. Algorithm 4.3 in [1].
!
! [1] Zhou, Y., Saad, Y., Tiago, M. L., Chelikowsky, J. R. (2006).
! Self-consistent-field calculations using Chebyshev-filtered subspace
! iteration. Journal of Computational Physics, 219(1), 172–184.
! http://doi.org/10.1016/j.jcp.2006.03.017
real(dp), intent(in) :: X(:,:)
real(dp), intent(out) :: Y(:,:)
procedure(mul_matrix_vector) :: HX ! returns y = H*x
integer, intent(in) :: m
real(dp), intent(in) :: a, b
real(dp), dimension(size(Y,1),size(Y,2)) :: X2, Ynew
real(dp) :: c, e, sigma, sigma1, sigma2
integer :: i, j
X2 = X
e = (b-a)/2; c = (b+a)/2 ! 1.
sigma = e/(a-c) ! 2.
sigma1 = sigma ! 3.
do j = 1, size(X2, 2)
    call HX(X2(:,j), Y(:,j)) ! 4.
end do
Y = (Y-c*X2)*sigma1/e ! 4.
do i = 2, m ! 5.
    sigma2 = 1/(2/sigma1-sigma) ! 6.
    do j = 1, size(X2, 2)
        call HX(X2(:,j), Ynew(:,j)) ! 7.
    end do
    Ynew = 2*(Ynew-c*Y)*sigma2/e - sigma*sigma2*X2 ! 7.
    X2 = Y ! 8.
    Y = Ynew ! 9.
    sigma = sigma2 ! 10.
end do ! 11.
end subroutine

subroutine solve_schroedinger2(myid, comm_all, commy, commz, Ng, nsub, Vloc, &
        L, G2, cutfn, nev, ncv, eigs, orbitals)
! Chebyshev filtering
integer, intent(in) :: myid, comm_all, commy, commz, Ng(3), nsub(3), nev, ncv
real(dp), intent(in) :: Vloc(:,:,:) ! Local effective potential
real(dp), intent(in) :: G2(:,:,:), L(:), cutfn(:,:,:)
real(dp), intent(out) :: eigs(:) ! eigs(nev)
! orbitals(Ng_local(1),Ng_local(2),Ng_local(3),nev)
real(dp), intent(inout) :: orbitals(:,:,:,:)

real(dp), allocatable :: d(:), v(:,:)
integer :: Ng_local(3), n
real(dp) :: t1, t2
logical :: verbose
verbose = .false.
Ng_local = Ng / nsub
n = product(Ng_local)
allocate(v(n,ncv), d(ncv))
call cpu_time(t1)
call peig(comm_all, myid, n, nev, ncv, "SA", av, d, v)
call cpu_time(t2)
eigs = d(:nev)+10
orbitals = reshape(v(:,:nev), [Ng_local(1),Ng_local(2),Ng_local(3),nev]) &
    * sqrt(product(Ng/L))
if (myid == 0 .and. verbose) then
    print *, "Arpack Time:", t2-t1
end if

contains

    subroutine av(x, y)
    ! Compute y = A*x
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    real(dp), dimension(Ng_local(1),Ng_local(2),Ng_local(3)) :: &
        psi
    psi = reshape(x, [Ng_local(1),Ng_local(2),Ng_local(3)])
    call applyH(commy, commz, Ng, nsub, Vloc-10, G2, cutfn, psi)
    y = reshape(psi, [product(Ng_local)])
    end

end subroutine


subroutine applyH(commy, commz, Ng, nsub, Vloc, G2, cutfn, psi)
! Apply the Hamiltonian: psi = H*psi
! The Hamiltonian is composed of the kinetic part and the local potential Vloc
! part. The `psi` and `Vloc` are given in real space.
integer, intent(in) :: commy, commz, Ng(:), nsub(:)
real(dp), intent(in) :: Vloc(:,:,:) ! Local potential in real space
real(dp), intent(in) :: G2(:,:,:), cutfn(:,:,:)
real(dp), intent(inout) :: psi(:,:,:) ! Wave function in real space
complex(dp), dimension(size(psi,1),size(psi,2),size(psi,3)) :: &
    psi2, psiG, psiG_vloc
call preal2fourier(psi, psiG, commy, commz, Ng, nsub)
psiG = psiG * cutfn
! psiG is our starting point

! Apply kinetic and potential
call pfourier2real(psiG, psi2, commy, commz, Ng, nsub)
call preal2fourier(Vloc*psi2, psiG_vloc, commy, commz, Ng, nsub)
psiG_vloc = psiG_vloc * cutfn

psiG = G2/2*psiG + psiG_vloc

! Convert to real space at the end
call pfourier2real(psiG, psi2, commy, commz, Ng, nsub)
psi = real(psi2,dp)
end subroutine


end module
