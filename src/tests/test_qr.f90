program test_qr
use types, only: dp
use linalg, only: qr_fact, eig, eigh
use utils, only: assert
implicit none
real(dp), allocatable :: A(:,:), Q(:,:), R(:,:)
real(dp), allocatable :: c(:,:), lam(:)
integer :: n, i
n = 3
allocate(A(n,n), Q(n,n), R(n,n), c(n,n), lam(n))
A(:,1) = [1, 2, 3]
A(:,2) = [2, 4, 2]
A(:,3) = [3, 2, 5]
call eigh(A, lam, c)
R = 0
do i = 1, n
    R(i,i) = lam(i)
end do
call assert(maxval(matmul(c, matmul(R, transpose(c))) - A) < 1e-12_dp)
print *, c(:,1)
call qr_fact(A, Q, R)
end
