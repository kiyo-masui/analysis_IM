module c_chol

use iso_c_binding, only: c_double, c_int, c_int64_t
use cholesky, only: cholesky

implicit none

contains

subroutine c_cholesky(n, a, lda) bind(c)
    integer(c_int64_t), intent(in):: n, lda
    real(c_double), dimension(lda,n), intent(inout) :: a
    integer :: nf, ldaf
    nf = n
    ldaf = lda
    call cholesky(n, a, lda)
end subroutine

end module
