subroutine cholesky(n,a,lda) !! a=l*l'; a=l
real, dimension(lda,n) :: a

call dpotrf('L',n,a,lda,info)
end
