      SUBROUTINE QCGETF(FOCK,N,FBLK,M,ISTART,IEND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FOCK(N,N),FBLK(M*M)
      IROW=0
      DO I=ISTART,IEND
      DO J=ISTART,IEND
         IROW=IROW+1
         FBLK(IROW)=FOCK(J,I)
      END DO
      END DO
      RETURN
      END 
