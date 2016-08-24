module choleskym

use iso_c_binding, only: c_double, c_int, c_int64_t


contains


!call test_cholesky(2000)
!end
! Cholesky and related subroutines.


subroutine c_cholesky(n, a, lda) bind(c)
    integer(c_int64_t), intent(in) :: n, lda
    real(c_double), dimension(lda,n), intent(inout) :: a
    integer :: nf, ldaf
    print *, n, lda
    nf = n
    ldaf = lda
    print *, nf, ldaf
    print *, a(1,1), a(1,2), a(2,2)
    call cholesky(nf, a(1,1), ldaf)
end subroutine

subroutine test_cholesky(n)
  integer n, i, t0, tr
  real, dimension(n,n) :: A
  real t_start, t_total
  
  do i=1,n
!    do j=i,n
!      A(j,i)=0.5
!    enddo
    A(i:,i) = 0.5*0
    A(i,i) = A(i,i) + i + 1
  enddo
  call system_clock(t0,tr)
  t_start = t0*1./tr
  call cholesky(n,A,n)
  write(*,*) a(:4,1)
  call system_clock(t0,tr)
  t_total = t0*1./tr
  t_total = t_total - t_start
  print*, 'Time was: ', t_total
end

subroutine tolower(a,n) !! set upper triangle to zeros
real a(n,n)
do i=1,n
        a(i,i+1:)=0
enddo
end

subroutine symm(a,n) !! copy lower triangle to upper triandgle
real a(n,n)
do i=1,n
        a(i,i+1:)=a(i+1:,i)
enddo
end


recursive subroutine cholesky(n,a,lda) !! a=l*l'; a=l
real, dimension(lda,n) :: a
! print*,n
 if (n .le. 1) then
     a(:n,:)=sqrt(a(:n,:))
     return
 endif
 n1=n/2; n2=n1+1; nd=n-n/2;
 call cholesky(n1,a(1,1),lda)
 call dtrsm(nd,n1,a(n2,1),lda,a(1,1),lda)
 call syrk(nd,n1,-1.,a(n2,1),lda,1.,a(n2,n2),lda)
 call cholesky(nd,a(n2,n2),lda)
end subroutine

recursive subroutine dtrsm(m,n,b,ldb,a,lda) ! b=b*(a^-1)'
real, dimension(ldb,n) :: b
real, dimension(lda,n) :: a
if (n .le. 1) then
    b(:m,:)=b(:m,:)*(1/a(1,1))
    return
endif
 n1=n/2; n2=n/2+1; nd=n-n/2
 call dtrsm(m,n1,b(1,1),ldb,a(1,1),lda)
 call dgemm('n','t',m,nd,n1,-1.,b(1,1),ldb,a(n2,1),lda,1.,b(1,n2),ldb)
 call dtrsm(m,nd,b(1,n2),ldb,a(n2,n2),lda)
end subroutine

recursive subroutine tinv(n,a,lda) ! a is lower triangle, a=a^-1
real, dimension(lda,n) :: a
if (n .le. 1) then
    a(:n,:)=1./a(:n,:)
    return
endif
 n1=n/2; n2=n/2+1; nd=n-n/2
 call tinv(n1,a,lda)
 call tinv(nd,a(n2,n2),lda)
 call trmm(nd,n1,a(n2,1),lda,a(1,1),lda)
 call trmml(nd,n1,a(n2,1),lda,a(n2,n2),lda)
end subroutine

recursive subroutine trmmt(m,n,b,ldb,a,lda) ! a is lower triangle, b=b*a'
real, dimension(ldb,n) :: b
real, dimension(lda,n) :: a
if (n .le. 1) then
    b(:m,:)=matmul(b(:m,:),transpose(a(:n,:)))
    return
endif
n1=n/2; n2=n/2+1; nd=n-n/2
 call trmmt(m,nd,b(1,n2),ldb,a(n2,n2),lda)
 call dgemm('n','t',m,nd,n1,1.,b,ldb,a(n2,1),lda,1.,b(1,n2),ldb)
 call trmmt(m,n1,b(1,1),ldb,a(1,1),lda)
end subroutine

recursive subroutine trmm(m,n,b,ldb,a,lda) ! a is lower triangle, b=-b*a !!!!!!!
real, dimension(ldb,n) :: b
real, dimension(lda,n) :: a
if (n .le. 1) then
    b(:m,:)=-matmul(b(:m,:),a(:n,:))      !!!!!!! NOTICE THE MINUS SIGN !!!!!!!
    return
endif
n1=n/2; n2=n/2+1; nd=n-n/2
 call trmm(m,n1,b(1,1),ldb,a(1,1),lda)
 call dgemm('n','n',m,n1,nd,-1.,b(1,n2),ldb,a(n2,1),lda,1.,b(1,1),ldb)
      !!!!!!! NOTICE THE MINUS SIGN !!!!!!!
 call trmm(m,nd,b(1,n2),ldb,a(n2,n2),lda)
end subroutine


recursive subroutine trmml(m,n,b,ldb,a,lda) ! a is lower triangle, b=a*b
real, dimension(ldb,n) :: b
real, dimension(lda,m) :: a
if (m .le. 1) then
    b(:m,:)=matmul(a(:m,:),b(:m,:))
    return
endif
m1=m/2; m2=m/2+1; md=m-m/2
 call trmml(md,n,b(m2,1),ldb,a(m2,m2),lda)
 call dgemm('n','n',md,n,m1,1.,a(m2,1),lda,b(1,1),ldb,1.,b(m2,1),ldb)
 call trmml(m1,n,b(1,1),ldb,a(1,1),lda)
end subroutine

recursive subroutine trmmlt(m,n,b,ldb,a,lda) ! a is lower triangle, b=a'*b
real, dimension(ldb,n) :: b
real, dimension(lda,m) :: a
 call tolower(a(:m,:m),m)
!b=matmul(transpose(a),b)
if (m .le. 1) then
    b(:m,:)=matmul(transpose(a(:m,:)),b(:m,:))
    return
endif
m1=m/2; m2=m/2+1; md=m-m/2
 call trmmlt(m1,n,b(1,1),ldb,a(1,1),lda)
 call dgemm('t','n',m1,n,md,1.,a(m2,1),lda,b(m2,1),ldb,1.,b(1,1),ldb)
 call trmmlt(md,n,b(m2,1),ldb,a(m2,m2),lda)
end subroutine

recursive subroutine slauum(n,a,lda) ! a is lower tri, a=a*a'
! returns lower triangle of c.
real, dimension(lda,n) :: a
if (n .le. 1) then
        a(:n,:)=matmul(a(:n,:),transpose(a(:n,:)))
        return
endif
n1=n/2; n2=n/2+1; nd=n-n/2
  call slauum(nd,a(n2,n2),lda)
  call syrk(nd,n1,1.,a(n2,1),lda,1.,a(n2,1),lda)
  call trmmt(nd,n1,a(n2,1),lda,a(1,1),lda)
  call slauum(n1,a(1,1),lda)
end

recursive subroutine slauumt(n,a,lda) ! a is lower tri, a=a'*a
! returns lower triangle of c.
real, dimension(lda,n) :: a
if (n .le. 1) then
        a(:n,:n)=matmul(transpose(a(:n,:n)),a(:n,:n))
        return
endif
n1=n/2; n2=n/2+1; nd=n-n/2
  call slauumt(n1,a(1,1),lda)
  call syrkl(n1,nd,1.,a(n2,1),lda,1.,a(1,1),lda)
  call trmmlt(nd,n1,a(n2,1),lda,a(n2,n2),lda)
  call slauumt(nd,a(n2,n2),lda)
end subroutine

recursive subroutine syrk(n,k,alpha,a,lda,beta,c,ldc) ! c=alpha a a'+beta c 
! returns lower triangle of c.
real, dimension(ldc,n) :: c
real, dimension(lda,k) :: a
if (n .le. 1) then
        c(:n,:)=alpha*matmul(a(:n,:),transpose(a(:n,:)))+beta*c(:n,:)
        return
endif
n1=n/2; n2=n/2+1; nd=n-n/2
  call syrk(n1,k,alpha,a(1,1),lda,beta,c(1,1),ldc)
  call syrk(nd,k,alpha,a(n2,1),lda,beta,c(n2,n2),ldc)
  call dgemm('n','t',nd,n1,k,alpha,a(n2,1),lda,a(1,1),lda,beta,c(n2,1),ldc)
end

recursive subroutine syrkl(n,k,alpha,a,lda,beta,c,ldc) ! c=alpha a' a+beta c 
! returns lower triangle of c.
real, dimension(ldc,n) :: c
real, dimension(lda,n) :: a                                       
if (n .le. 1) then
        c(:k,:)=alpha*matmul(transpose(a(:k,:)),a(:k,:))+beta*c(:k,:)
        return
endif
n1=n/2; n2=n/2+1; nd=n-n/2
  call syrkl(n1,k,alpha,a(1,1),lda,beta,c(1,1),ldc)
  call syrkl(nd,k,alpha,a(1,n2),lda,beta,c(n2,n2),ldc)
  call dgemm('t','n',nd,n1,k,alpha,a(1,n2),lda,a(1,1),lda,beta,c(n2,1),ldc)
end

end module
