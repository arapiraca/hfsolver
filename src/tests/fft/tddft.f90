program tddft
use types, only: dp
use constants, only: Ha2eV, density2gcm3, u2au, s2au, K2au, i_, pi, &
    ang2bohr, bohr2ang
use fourier, only: dft, idft, fft, fft_vectorized, fft_pass, fft_pass_inplace, &
        fft_vectorized_inplace, calculate_factors, ifft_pass, fft2_inplace, &
        fft3_inplace, ifft3_inplace
use utils, only: assert, init_random, stop_error, get_int_arg, get_float_arg, &
    allocate_mold, clock, whitechar
use ffte, only: factor
use ofdft, only: read_pseudo
use pofdft_fft, only: pfft3_init, preal2fourier, pfourier2real, &
    real_space_vectors, reciprocal_space_vectors, calculate_myxyz, &
    pintegral, pintegralG, free_energy, free_energy_min, &
    radial_potential_fourier, psum, pmaxval, distribute, poisson_kernel, &
    collate
use openmp, only: omp_get_wtime
use mpi2, only: mpi_finalize, MPI_COMM_WORLD, mpi_comm_rank, &
    mpi_comm_size, mpi_init, mpi_comm_split, MPI_INTEGER, &
    mpi_barrier, mpi_bcast, MPI_DOUBLE_PRECISION, MPI_LOGICAL
use md, only: positions_bcc, positions_fcc
use arpack, only: peig, eig
use pksdft_fft, only: solve_schroedinger, applyH
use xc, only: xc_pz
use mixings, only: mixing_linear, mixing_anderson
use ewald_sums, only: ewald_box
use efermi, only: fermi_dirac_smearing
implicit none

complex(dp), dimension(:,:,:), allocatable :: neG, psiG, psi
real(dp), dimension(:,:,:), allocatable :: G2, Htot, HtotG, Ven0G, ne, Vloc, &
    Veff
real(dp), allocatable :: G(:,:,:,:), X(:,:,:,:), Xion(:,:), q(:), &
    current(:,:,:,:), eigs(:), orbitals(:,:,:,:), occ(:), &
    Vee(:,:,:), Vee_xc(:,:,:), exc(:,:,:), Vxc(:,:,:), forces(:,:), &
    cutfn(:,:,:), cutfn2(:,:,:), tmp_global(:,:,:), psi_r(:,:,:)
complex(dp), allocatable :: dpsi(:,:,:,:), VeeG(:,:,:), VenG(:,:,:), &
    corbitals(:,:,:,:), tmp(:,:,:)
real(dp) :: L(3), stress(6), current_avg(3), A
integer :: i, j, k, u, u2, u3
integer :: Ng(3)
integer :: natom, nelec, nband, max_iter, field_dir
logical :: velocity_gauge
real(dp) :: T_au, dt, rho, norm, w2, Vmin, Ekin, Etot, &
    Eee, Een_loc, E_xc, Enn, Een_core, G2cut, G2cut2
real(dp) :: rloc, C1, C2, Zion, Ecut, E0, Ex, t, td, tw
real(dp), allocatable :: m(:)
integer :: it

!  parallel variables
integer :: comm_all, commy, commz, nproc, ierr, nsub(3), Ng_local(3)
integer :: myid ! my ID (MPI rank), starts from 0
integer :: myxyz(3) ! myid, converted to the (x, y, z) box, starts from 0


call mpi_init(ierr)
comm_all  = MPI_COMM_WORLD
call mpi_comm_rank(comm_all, myid, ierr)
call mpi_comm_size(comm_all, nproc, ierr)

if (myid == 0) then
    call load_initial_pos(natom, L, Xion)
end if
call mpi_bcast(natom, 1, MPI_INTEGER, 0, comm_all, ierr)
if (myid /= 0) then
    allocate(Xion(3,natom))
end if
call mpi_bcast(Xion, size(Xion), MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
call mpi_bcast(L, size(L), MPI_DOUBLE_PRECISION, 0, comm_all, ierr)

if (myid == 0) then
    call read_params(Ng, T_au, dt, Ecut, nband)
    call read_input(nsub, dt, max_iter, E0, td, tw, velocity_gauge, field_dir)
end if
call mpi_bcast(Ng, size(Ng), MPI_INTEGER, 0, comm_all, ierr)
call mpi_bcast(nband, 1, MPI_INTEGER, 0, comm_all, ierr)
call mpi_bcast(max_iter, 1, MPI_INTEGER, 0, comm_all, ierr)
call mpi_bcast(nsub, size(nsub), MPI_INTEGER, 0, comm_all, ierr)
call mpi_bcast(T_au, 1, MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
call mpi_bcast(dt, 1, MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
call mpi_bcast(Ecut, 1, MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
call mpi_bcast(E0, 1, MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
call mpi_bcast(td, 1, MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
call mpi_bcast(tw, 1, MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
call mpi_bcast(velocity_gauge, 1, MPI_LOGICAL, 0, comm_all, ierr)
call mpi_bcast(field_dir, 1, MPI_INTEGER, 0, comm_all, ierr)

G2cut = 2*Ecut
G2cut2 = 4*G2cut
call assert(all(Ng /= -1))

Ng_local = Ng / nsub


allocate(m(natom))
nelec = natom
m = 2._dp * u2au ! Using D mass in atomic mass units [u]

!L = (sum(m) / rho)**(1._dp/3)
rho = sum(m) / product(L)

allocate(q(natom))
q = 1

if (myid == 0) then
    print *, "Input:"
    print *, "N =", natom
    print *, "L =", L, "a.u. =", L * bohr2ang, "Angst"
    print *, "T =", T_au / K2au, "K =", T_au, "a.u. =", T_au * Ha2eV, "eV"
    print "('dt =', es10.2, ' s = ', es10.2, ' fs = ', es10.2, ' a.u.')", &
        dt/s2au, dt/s2au * 1e15_dp, dt
    print *, "Ecut =", Ecut, "a.u. =", Ecut * Ha2eV, "eV"
    print *, "nband =", nband
    print *, "E0 =", E0
    print *, "td =", td
    print *, "tw =", tw
    print *, "velocity_gauge =", velocity_gauge
    print *
    print *, "Calculated quantities:"
    print *, "rho = ", rho * density2gcm3, "g/cc = ", rho, "a.u."
    print *
    print *, "nproc:   ", nproc
    print *, "nsub:    ", nsub
    print *, "Ng:      ", Ng
    print *, "Ng_local:", Ng_local

    if (product(nsub) /= nproc) then
        call stop_error("nproc must be equal to the number of subdomains")
    end if

end if

call pfft3_init(myid, comm_all, Ng, nsub)

myxyz = calculate_myxyz(myid, nsub)

! Note that myxyz(3) corresponds to commy, and myxyz(2) to commz
call mpi_comm_split(comm_all, myxyz(3), 0, commy, ierr)
call mpi_comm_split(comm_all, myxyz(2), 0, commz, ierr)


! We allocate the global array here, to ensure we have enough memory
allocate(tmp_global(Ng(1), Ng(2), Ng(3)))

allocate(ne(Ng_local(1), Ng_local(2), Ng_local(3)))
allocate(neG(Ng_local(1), Ng_local(2), Ng_local(3)))
call allocate_mold(G2, ne)
call allocate_mold(Htot, ne)
call allocate_mold(HtotG, ne)
call allocate_mold(Vloc, ne)
call allocate_mold(Veff, ne)
call allocate_mold(Vxc, ne)
call allocate_mold(exc, ne)
call allocate_mold(Ven0G, ne)
call allocate_mold(Vee, ne)
call allocate_mold(Vee_xc, ne)
call allocate_mold(psi_r, ne)
call allocate_mold(psiG, neG)
call allocate_mold(psi, neG)
call allocate_mold(VeeG, neG)
call allocate_mold(VenG, neG)
allocate(orbitals(Ng_local(1), Ng_local(2), Ng_local(3), nband))
allocate(corbitals(Ng_local(1), Ng_local(2), Ng_local(3), nband))
allocate(X(Ng_local(1), Ng_local(2), Ng_local(3), 3))
allocate(dpsi(Ng_local(1), Ng_local(2), Ng_local(3), 3))
call allocate_mold(tmp, psiG)
call allocate_mold(G, X)
call allocate_mold(current, X)
allocate(forces(3, natom))
! For now assume a box, until positions_bcc can accept a vector L(:)
! And radial_potential_fourier
call assert(abs(L(2)-L(1)) < 1e-15_dp)
call assert(abs(L(3)-L(1)) < 1e-15_dp)
!call positions_fcc(Xion, L(1))
!Xion(:,1) = [-0.7_dp+L(1)/2, L(2)/2, L(3)/2]
!Xion(:,2) = [-0.7_dp+L(1)/2, L(2)/2, L(3)/2]
!Xion(:,1) = L/2
!Xion = 0
!Xion(1,1) = L(1)/2
!Xion = Xion + 1e-4_dp ! Shift the atoms, so that 1/R is defined
call real_space_vectors(L, X, Ng, myxyz)
call reciprocal_space_vectors(L, G, G2, Ng, myxyz)

call ewald_box(L(1), Xion, q, Enn, forces, stress)

! HDH pseudopotential for Hydrogen
rloc = 0.2_dp
C1 = -4.180237_dp
C2 =  0.725075_dp
Zion = 1
! HDH local pseudopotential for Pb
!rloc = 0.617500_dp
!C1 =   0.753143_dp
!C2 =   0
!Zion = 4

do k = 1, Ng_local(3)
do j = 1, Ng_local(2)
do i = 1, Ng_local(1)
    w2 = G2(i,j,k)
    if (w2 < tiny(1._dp)) then
        Ven0G(i,j,k) = 0
    else
        Ven0G(i,j,k) = -4*pi*Zion/w2*exp(-w2*rloc**2/2) &
            + sqrt(8*pi**3)*rloc**3*exp(-w2*rloc**2/2)*(C1+C2*(3-w2*rloc**2))
    end if
end do
end do
end do
Ven0G = Ven0G / product(L)


VenG = 0
do i = 1, natom
    VenG = VenG + Ven0G * exp(-i_ * &
        (G(:,:,:,1)*Xion(1,i) + G(:,:,:,2)*Xion(2,i) + G(:,:,:,3)*Xion(3,i)))
end do

!Vloc = 0
!do na = 1, natom
!    do k = 1, Ng_local(3)
!    do j = 1, Ng_local(2)
!    do i = 1, Ng_local(1)
!        r = sqrt(sum((X(i,j,k,:)-Xion(:,na))**2))
!        Vloc(i,j,k) = Vloc(i,j,k) - Zion/r * erf(r/(sqrt(2._dp)*rloc)) &
!            + exp(-1._dp/2*(r/rloc)**2) * (C1 + C2*(r/rloc)**2)
!    end do
!    end do
!    end do
!end do

Vmin = minval(Vloc)

!if (myid == 0) then
!    open(newunit=u, file="a.txt", status="replace")
!    write(u,*) Vloc(:,16,16)
!end if

!print *, "Ecut:", Ecut
!print *, "G2cut:", G2cut / (2*pi)**2
!print *, "boxsq:", maxval(G2(:,1,1)) / (2*pi)**2
!print *, "boxcut:", sqrt(maxval(G2(:,1,1)) / G2cut)
!print *, "gsqcut:", G2cut2 / (2*pi)**2

allocate(cutfn(Ng_local(1),Ng_local(2),Ng_local(3)))
allocate(cutfn2(Ng_local(1),Ng_local(2),Ng_local(3)))
where (G2 > G2cut)
    cutfn = 0
elsewhere
    !cutfn = (1-G2/G2cut)**12
    cutfn = 1
end where
!print *, maxval(G2), G2cut, 2.04204_dp**2 * G2cut
!stop "OK"
where (G2 > G2cut2)
    cutfn2 = 0
elsewhere
    cutfn2 = 1
end where

!call pfourier2real(G2*VenG/(4*pi), Vloc, commy, commz, Ng, nsub)
VenG = VenG * cutfn2
call pfourier2real(VenG, Vloc, commy, commz, Ng, nsub)
norm = pintegral(comm_all, L, Vloc, Ng)
if (myid == 0) then
    print *, "INT:", norm
end if

if (myid == 0) print *, "Loading orbitals"
if (myid == 0) open(newunit=u, file="orbitals.dat", status="old")
do i = 1, nband
    if (myid == 0) read(u,*) tmp_global
    call distribute(comm_all, myid, nsub, 0, tmp_global, orbitals(:,:,:,i))
end do
if (myid == 0) close(u)
if (myid == 0) print *, "Square of norms of orbitals <psi|psi>:"
do i = 1, nband
    norm = pintegral(comm_all, L, orbitals(:,:,:,i)**2, Ng)
    if (myid == 0) print *, i, norm
end do

if (myid == 0) print *, "Loading Veff"
if (myid == 0) open(newunit=u, file="Veff.dat", status="old")
if (myid == 0) read(u,*) tmp_global
call distribute(comm_all, myid, nsub, 0, tmp_global, Veff)
if (myid == 0) close(u)

if (myid == 0) print *, "Compute eigs using a Hamiltonian and compare " // &
    "against KS reference results:"
if (myid == 0) open(newunit=u, file="eigs.dat", status="old")
allocate(eigs(nband), occ(nband))
if (myid == 0) read(u,*) eigs
if (myid == 0) read(u,*) occ
if (myid == 0) close(u)
call mpi_bcast(eigs, size(eigs), MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
call mpi_bcast(occ, size(occ), MPI_DOUBLE_PRECISION, 0, comm_all, ierr)

if (myid == 0) print "(a6, a16, a10, a10)", "n", "<psi|H|psi>", "Error", "occ"
do i = 1, nband
    psi_r = orbitals(:,:,:,i)
    call applyH(commy, commz, Ng, nsub, Veff, G2, cutfn, psi_r)
    norm = pintegral(comm_all, L, orbitals(:,:,:,i)*psi_r, Ng)
    if (myid == 0) print "(i6, f16.10, es10.2, f10.5)", &
        i, norm, abs(norm - eigs(i)), occ(i)
end do

if (myid == 0) print *, "Propagation"

if (myid == 0) open(newunit=u, file="energies.txt", status="replace")
if (myid == 0) open(newunit=u2, file="current.txt", status="replace")
if (myid == 0) open(newunit=u3, file="current_density.txt", status="replace")
corbitals = orbitals
t = 0
do it = 1, max_iter
    t = t + dt
    if (myid == 0) print *, "Starting Iteration:", it

    Ex = E0 * exp(-(t-td)**2/(2*tw**2)) / (sqrt(2*pi)*tw)
    A = -E0 * (erf(sqrt(2._dp)*(t - td)/(2*tw))/2 + 1._dp/2)
    call preal2fourier(Veff, psiG, commy, commz, Ng, nsub)
    psiG = psiG * cutfn2
    call pfourier2real(psiG, Veff, commy, commz, Ng, nsub)
    if (velocity_gauge) then
        ! velocity gauge
        Htot = Veff+A**2/2
        HtotG = G2/2 + A*G(:,:,:,field_dir)
    else
        ! length gauge
        Htot = Veff+X(:,:,:,field_dir)*Ex
        HtotG = G2/2
    end if
    do i = 1, nband
        psi = corbitals(:,:,:,i) * exp(-i_*Htot*dt/2)
        call preal2fourier(psi, psiG, commy, commz, Ng, nsub)
        psiG = psiG * cutfn
        psiG = psiG * exp(-i_*HtotG*dt)
        call pfourier2real(psiG, psi, commy, commz, Ng, nsub)
        psi = psi * exp(-i_*Htot*dt/2)

        ! Apply a cutoff
        call preal2fourier(psi, psiG, commy, commz, Ng, nsub)
        psiG = psiG * cutfn
        call pfourier2real(psiG, psi, commy, commz, Ng, nsub)

        ! Normalize
        norm = pintegral(comm_all, L, abs(psi)**2, Ng)
        psi = psi / sqrt(norm)

        corbitals(:,:,:,i) = psi
    end do
    if (myid == 0) print *, "Square of norms of orbitals <psi|psi>:"
    do i = 1, nband
        norm = pintegral(comm_all, L, abs(corbitals(:,:,:,i))**2, Ng)
        if (myid == 0) print *, i, norm
    end do

    ! Density
    ne = 0
    Ekin = 0
    do i = 1, nband
        ne = ne + occ(i)*abs(corbitals(:,:,:,i))**2
        call preal2fourier(corbitals(:,:,:,i), psiG, commy, commz, Ng, nsub)
        Ekin = Ekin &
            + occ(i) * 1._dp/2 * pintegralG(comm_all, L, G2*abs(psiG)**2)
    end do

    ! Hartree (Poisson)
    call preal2fourier(ne, neG, commy, commz, Ng, nsub)
    call poisson_kernel(myid, size(neG), neG, G2, VeeG)
    VeeG = VeeG * cutfn2
    call pfourier2real(VeeG, Vee, commy, commz, Ng, nsub)

    ! XC
    call xc_pz(ne, exc, Vxc)

    Een_loc = pintegral(comm_all, L, Vloc*ne, Ng)
    E_xc = pintegral(comm_all, L, exc*ne, Ng)
    Eee = pintegral(comm_all, L, Vee*ne, Ng) / 2
    Een_core = 0 ! For now

    Veff = Vloc + Vee + Vxc
    Etot = Ekin + Eee + E_xc + Enn + Een_core + Een_loc

    if (myid == 0) then
        print *, "n, E[a.u.] E[eV] occ"
        print "(a, es22.14)", "Ekin:     ", Ekin
        print "(a, es22.14)", "Eee:      ", Eee
        print "(a, es22.14)", "Exc:      ", E_xc
        print "(a, es22.14)", "Enn:      ", Enn
        print "(a, es22.14)", "Een_core: ", Een_core
        print "(a, es22.14)", "Een_loc:  ", Een_loc
        print "(a, es22.14)", "Een_NL:   ", 0._dp
        print "(a, es22.14)", "Etot:     ", Etot

        write(u,*) t, Etot, Ekin, Eee, E_xc, Enn, Een_core, Een_loc, &
            Ex, A
    end if

    ! Current
    current = 0
    do i = 1, nband
        psi = corbitals(:,:,:,i)
        call preal2fourier(corbitals(:,:,:,i), psiG, commy, commz, Ng, nsub)
        do j = 1, 3
            call pfourier2real(i_*G(:,:,:,j)*psiG, dpsi(:,:,:,j), &
                    commy, commz, Ng, nsub)
            tmp = i_/(2*natom) * (conjg(psi)*dpsi(:,:,:,j)- &
                psi*conjg(dpsi(:,:,:,j)))
            if (velocity_gauge) then
                if (j == field_dir) tmp = tmp - A*ne/natom
            end if
            if (maxval(abs(aimag(tmp))) > 1e-12_dp) then
                print *, "INFO: current  max imaginary part:", &
                    maxval(aimag(tmp))
            end if
            current(:,:,:,j) = current(:,:,:,j) + occ(i) * real(tmp,dp)
        end do
    end do
    do j = 1, 3
        current_avg(j) = pintegral(comm_all, L, current(:,:,:,j),Ng)/product(L)
    end do
    if (myid == 0) write(u2,*) t, current_avg
    if (mod(it, 10) == 0) then
        call collate(comm_all, myid, nsub, 0, current(:,:,:,1), tmp_global)
        if (myid == 0) write(u3,*) tmp_global(:,:,Ng(3)/2)
        call collate(comm_all, myid, nsub, 0, current(:,:,:,2), tmp_global)
        if (myid == 0) write(u3,*) tmp_global(:,:,Ng(3)/2)
    end if
end do

if (myid == 0) close(u)
if (myid == 0) close(u2)
if (myid == 0) close(u3)
if (myid == 0) print *, "Done"

call mpi_finalize(ierr)

contains

    subroutine read_input(nsub, dt, max_iter, E0, td, tw, velocity_gauge, &
            field_dir)
    integer, intent(out) :: nsub(3), max_iter, field_dir
    real(dp), intent(out) :: dt ! in a.u.
    real(dp), intent(out) :: E0, td, tw
    logical, intent(out) :: velocity_gauge
    integer :: LNPU(3)
    namelist /domain/ nsub, dt, max_iter, E0, td, tw, velocity_gauge, field_dir
    integer :: u
    nsub = -1
    dt = -1
    max_iter = -1
    velocity_gauge = .true.
    E0 = -1
    td = -1
    tw = -1
    field_dir = 1
    open(newunit=u, file="input_tddft", status="old")
    read(u,nml=domain)
    close(u)
    if (nsub(1) == -1) then
        call factor(nproc, LNPU)
        nsub(3) = (2**(LNPU(1)/2))*(3**(LNPU(2)/2))*(5**(LNPU(3)/2))
        nsub(2) = nproc / nsub(3)
        nsub(1) = 1
    end if
    if (dt < 0) call stop_error("dt is not specified")
    if (max_iter < 0) call stop_error("max_iter is not specified")
    if (E0 < 0) call stop_error("E0 not specified")
    if (td < 0) call stop_error("td not specified")
    if (tw < 0) call stop_error("tw not specified")
    end subroutine

    subroutine read_params(Ng, T, dt, Ecut, nband)
    integer, intent(out) :: Ng(:), nband
    real(dp), intent(out) :: T, dt, Ecut  ! in a.u.
    namelist /domain/ Ng, T, dt, Ecut, nband
    integer :: u
    open(newunit=u, file="params.out", status="old")
    read(u,nml=domain)
    close(u)
    end subroutine

    subroutine load_initial_pos(natom, L, Xion)
    integer, intent(out) :: natom
    real(dp), intent(out) :: L(:)
    real(dp), intent(out), allocatable :: Xion(:,:)
    real(dp) :: A, t(3,3)
    integer :: u, natom_types
    character :: c
    integer, allocatable :: ncounts(:)
    open(newunit=u, file="POSCAR", status="old")
    read(u,*) ! 1 line description
    read(u,*) A
    read(u,*) t(:, 1)
    read(u,*) t(:, 2)
    read(u,*) t(:, 3)

    t = t*A*ang2bohr
    L(1) = t(1,1)
    L(2) = t(2,2)
    L(3) = t(3,3)
    call assert(abs(t(2, 1)) < tiny(1._dp))
    call assert(abs(t(3, 1)) < tiny(1._dp))
    call assert(abs(t(1, 2)) < tiny(1._dp))
    call assert(abs(t(3, 2)) < tiny(1._dp))
    call assert(abs(t(1, 3)) < tiny(1._dp))
    call assert(abs(t(2, 3)) < tiny(1._dp))

    natom_types = number_of_columns(u)
    allocate(ncounts(natom_types))
    read(u,*) ! Atom types
    read(u,*) ncounts ! Atom numbers per type
    natom = sum(ncounts)
    allocate(Xion(3,natom))
    read(u,*) c ! Selective dynamics
    call assert(c == "S")
    read(u,*) c ! Direct
    call assert(c == "D")
    do i = 1, natom
        read(u,*) Xion(:,i)
        Xion(:,i) = Xion(:,i)*L
    end do
    close(u)
    end subroutine

    integer function number_of_columns(u) result(ncol)
    ! Determines the number of columns on a line (does not affect file
    ! position).
    integer, intent(in) :: u ! file handle
    integer :: ios
    character :: c
    logical :: lastwhite
    ncol = 0
    lastwhite = .true.
    do
       read(u, '(a)', advance='no', iostat=ios) c
       if (ios /= 0) exit
       if (lastwhite .and. .not. whitechar(c)) ncol = ncol + 1
       lastwhite = whitechar(c)
    end do
    backspace(u)
    end function

end program
