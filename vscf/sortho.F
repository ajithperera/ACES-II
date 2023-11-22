      SUBROUTINE SORTHO(C,SOVRLP,SCR1,SCR2,SCR3,LDIM2,NBAS)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION C(NBAS*NBAS),SOVRLP(NBAS*NBAS)
      DIMENSION SCR1(LDIM2),SCR2(LDIM2),SCR3(LDIM2)
C
      DIMENSION JFLAG(16)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
#include "symm2.com"
      COMMON /FLAGS/ IFLAGS(100)
C
      DATA ZILCH /0.0/
      DATA ONEM /-1.0/
      DATA ONE /1.0/
      DATA TWO /2.0/
      DATA I16 /16/
C
      INDX2(I,J,N)=I+(J-1)*N
C
      WRITE(LUOUT,1000)
 1000 FORMAT('  @SORTHO-I, Orthonormalizing initial guess. ')
C
C     Subroutine to orthonormailize a set of eigenvectors.
C
C     Strategy :
C                  T
C                 C  S C = X
C
C     C is not orthonormal set, but C * X **(-1/2) is. Hence we
C     diagonalize X to form X**(-1/2).
C
C     Read overlap matrix.
C
      CALL GETREC(20,'JOBARC','AOOVRLAP',NBAS*NBAS*IINTFP,SOVRLP)

      DO   100 I=1,NIRREP
      IF(NBFIRR(I).EQ.0) GOTO 100
C
      CALL GETBLK(SOVRLP,SCR1,NBFIRR(I),NBAS,IREPS(I))
      CALL GETBLK(C     ,SCR2,NBFIRR(I),NBAS,IREPS(I))
C
      CALL ZERO(SCR3,NBFIRR(I)*NBFIRR(I))
      CALL XGEMM('N','N',NBFIRR(I),NBFIRR(I),NBFIRR(I), 1.0D+00,
     1           SCR1,NBFIRR(I),SCR2,NBFIRR(I), 0.0D+00,
     1           SCR3,NBFIRR(I))
      CALL ZERO(SCR1,NBFIRR(I)*NBFIRR(I))
      CALL XGEMM('T','N',NBFIRR(I),NBFIRR(I),NBFIRR(I), 1.0D+00,
     1           SCR2,NBFIRR(I),SCR3,NBFIRR(I), 0.0D+00,
     1           SCR1,NBFIRR(I))
      
C
C     Ith symmetry block of X is in SCR1. Diagonalize :
C
C     Eigenvectors in SCR2, eigenvalues in SCR1.
C
      CALL ZERO(SCR2,NBFIRR(I)*NBFIRR(I))
      CALL  EIG(SCR1,SCR2,NBFIRR(I),NBFIRR(I),0)
C
C     Raise diagonal matrix to power of -1/2.
C
      DO   20 K=1,NBFIRR(I)
      DO   10 J=1,NBFIRR(I)
      IF(J.NE.K )THEN
      SCR1(INDX2(J,K,NBFIRR(I))) = 0.0D+00
      ELSE
          SCR1(INDX2(J,J,NBFIRR(I))) = 1.0D+00 / 
     1                           DSQRT(SCR1(INDX2(J,J,NBFIRR(I))))
      ENDIF
   10 CONTINUE
   20 CONTINUE
C
C     Compute X**-1/2 = U * DIAG**-1/2 * U(T).
C
      CALL XGEMM('N','T',NBFIRR(I),NBFIRR(I),NBFIRR(I),1.0D+00,
     1           SCR1,NBFIRR(I),SCR2,NBFIRR(I),0.0D+00,
     1           SCR3,NBFIRR(I))
      CALL XGEMM('N','N',NBFIRR(I),NBFIRR(I),NBFIRR(I),1.0D+00,
     1           SCR2,NBFIRR(I),SCR3,NBFIRR(I),0.0D+00,
     1           SCR1,NBFIRR(I))
C
C     Replace original MO coefficients by orthonormal set, C*X**-1/2.
C
      CALL GETBLK(C     ,SCR2,NBFIRR(I),NBAS,IREPS(I))
      CALL  XGEMM('N','N',NBFIRR(I),NBFIRR(I),NBFIRR(I),1.0D+00,
     1            SCR2,NBFIRR(I),SCR1,NBFIRR(I),0.0D+00,
     1            SCR3,NBFIRR(I))
C
      CALL PUTBLK(C,SCR3,NBFIRR(I),NBAS,IREPS(I))
C
C     Compute X now with new set. C is in SCR3.
C
      CALL GETBLK(SOVRLP,SCR2,NBFIRR(I),NBAS,IREPS(I))
      CALL XGEMM('N','N',NBFIRR(I),NBFIRR(I),NBFIRR(I), 1.0D+00,
     1           SCR2,NBFIRR(I),SCR3,NBFIRR(I), 0.0D+00,
     1           SCR1,NBFIRR(I))
      CALL XGEMM('T','N',NBFIRR(I),NBFIRR(I),NBFIRR(I), 1.0D+00,
     1           SCR3,NBFIRR(I),SCR1,NBFIRR(I), 0.0D+00,
     1           SCR2,NBFIRR(I))
C
C      IF(IFLAGS(1).GE.5)THEN
C      CALL OUTPUT(SCR2,1,NBFIRR(I),1,NBFIRR(I),NBFIRR(I),NBFIRR(I),1)
C      ENDIF
C                 
  100 CONTINUE
      RETURN
      END
