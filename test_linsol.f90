program test_linsol
implicit none
integer, parameter::n=5
real::A(n,n),b(n),x(n)
A = reshape(  (/1., 2.,0., 0.,0.,&
                2., 3.,1., 0.,0.,&
                0.,-3.,4., 2.,0.,&
                0., 0.,4., 7.,1.,&
                0., 0.,0.,-5.,6./),&
                (/5,5/) )
A = transpose(A)
b = (/5,9,2,19,-4/)
! print*, transpose(A)
! print*, b
call linsol(A,b,x,n)
! print*, transpose(A)
! print*, b
! print*,
print*, x ! correct solution is : 1,2,1,2,1
endprogram

subroutine linsol(A, b, x, nnod)
    implicit none

    integer, intent(in):: nnod
    real, intent(inout)::A(nnod,nnod), b(nnod)
    real, intent(inout)::x(nnod)
    integer :: i 
    real :: w


    do i = 2, nnod
        w = A(i,i-1)/A(i-1,i-1)
        ! print*, w, A(i,i-1), A(i-1,i-1)
        A(i,i) = A(i,i) - w * A(i-1,i)
        b(i) = b(i) - w * b(i-1)
        ! print*,i,'w', w, 'A', A(i,i), 'b', b(i)
        ! print*, transpose(A)
    enddo
    x(nnod) = b(nnod) / A(nnod, nnod)
    do i = nnod-1, 1, -1
        x(i) = (b(i) - A(i,i+1) * x(i+1) ) / A(i,i)
    enddo
endsubroutine

SUBROUTINE SMLINN(A,X,B,NMX,N)
    IMPLICIT NONE
    INTEGER NMX,N
    REAL A(NMX,NMX),X(NMX),B(NMX)
    REAL R
    INTEGER K,I,J
! C     Form X = A^{-1} B
! C     Useful subroutine for inverse
! C     This sub overwrites the matrix A. 
    DO K=1,N-1
       DO I=K+1,N
          A(I,K)=A(I,K)/A(K,K)
       END DO
       DO J=K+1,N
          DO I=K+1,N
             A(I,J)=A(I,J) - A(I,K)*A(K,J)
          END DO
       END DO
    END DO
! C     
! C     Solve L_1 x=b
    DO I=1,N
       R=0.
       DO J=1,I-1
          R=R+A(I,J)*X(J)
       END DO
       X(I)=B(I)-R
    END DO
! C     
! C     Solve U x=y
    DO I=N,1,-1
       R=0.
       DO J=I+1,N
          R=R+A(I,J)*X(J)
       END DO
       X(I)=(X(I)-R)/A(I,I)
    END DO
    RETURN
    END