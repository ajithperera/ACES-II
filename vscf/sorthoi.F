      SUBROUTINE SORTHOI(C,SOVRLP,SCR1,SCR2,SCR3,LDIM2,NBAS,ALLMOS,NOCC,
     &                   IRREP)
C
C     ALLMOS --- Indicates whether we are orthonormalizing all MOs
C                or just the occupied MOs. In either case, C has the
C                dimensions and ordering of the full set of MOs (ie it
C                is an NBAS by NBAS matrix). If ALLMOS is .FALSE., it is
C                best if the only nonzero elements of C are the occupied
C                orbitals.
C     NOCC   --- Contains occupation for current spin.
C     NMOS   --- The number of orbitals to be orthogonalized in the
C                current symmetry block.
C     IRREP  --- Doing just one irrep at once.
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C-----------------------------------------------------------------------
      DOUBLE PRECISION C,SOVRLP,SCR1,SCR2,SCR3
      INTEGER LDIM2,NBAS,NOCC
      LOGICAL ALLMOS
C-----------------------------------------------------------------------
      DIMENSION C(NBAS*NBAS),SOVRLP(NBAS*NBAS)
      DIMENSION SCR1(LDIM2),SCR2(LDIM2),SCR3(LDIM2),NOCC(8)
C-----------------------------------------------------------------------
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/  LUOUT,MOINTS
#include "symm2.com"
      COMMON /FLAGS/  IFLAGS(100)
C
      DATA ZILCH /0.0/
      DATA ONEM /-1.0/
      DATA ONE /1.0/
C
      INDX2(I,J,N)=I+(J-1)*N
C
      WRITE(LUOUT,1000)
 1000 FORMAT(' @SORTHOI-I, Orthonormalizing initial guess. ')
C
C     Subroutine to orthonormalize a set of eigenvectors.
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
C
ccc      DO   100 I=1,NIRREP
C
      IF(ALLMOS)THEN
       NMOS = NBFIRR(IRREP)
      ELSE
       NMOS = NOCC(IRREP)
      ENDIF
C
ccc      IF(NMOS.EQ.0) GOTO 100
C
      CALL ZERO(SCR1,LDIM2)
      CALL ZERO(SCR2,LDIM2)
      CALL ZERO(SCR3,LDIM2)
C
      CALL GETBLK(SOVRLP,SCR1,NBFIRR(IRREP),NBAS,IREPS(IRREP))
ccc      CALL GETBLK(C     ,SCR2,NBFIRR(I),NBAS,IREPS(I))
c     CALL ICOPY(C,SCR2,NBFIRR(IRREP)*NBFIRR(IRREP)*IINTFP)
      CALL ICOPY(NBFIRR(IRREP)*NBFIRR(IRREP)*IINTFP,C,1,SCR2,1)
C
      CALL XGEMM('N','N',NBFIRR(IRREP),NMOS,NBFIRR(IRREP), 1.0D+00,
     &           SCR1,NBFIRR(IRREP),SCR2,NBFIRR(IRREP), 0.0D+00,
     &           SCR3,NBFIRR(IRREP))
      CALL ZERO(SCR1,LDIM2)
      CALL XGEMM('T','N',NMOS,NMOS,NBFIRR(IRREP), 1.0D+00,
     &           SCR2,NBFIRR(IRREP),SCR3,NBFIRR(IRREP), 0.0D+00,
     &           SCR1,NMOS)
C
C     Ith symmetry block of X is in SCR1. Diagonalize :
C
C     Eigenvectors in SCR2, eigenvalues in SCR1.
C
      CALL ZERO(SCR2,LDIM2)
      CALL  EIG(SCR1,SCR2,NMOS,NMOS,0)
C
C     Raise diagonal matrix to power of -1/2.
C
      DO   20 K=1,NMOS
      DO   10 J=1,NMOS
      IF(J.NE.K)THEN
       SCR1(INDX2(J,K,NMOS)) = 0.0D+00
      ELSE
       SCR1(INDX2(J,J,NMOS)) = 1.0D+00 / 
     &                         DSQRT(SCR1(INDX2(J,J,NMOS)))
      ENDIF
   10 CONTINUE
   20 CONTINUE
C
C     Compute X**-1/2 = U * DIAG**-1/2 * U(T).
C
      CALL XGEMM('N','T',NMOS,NMOS,NMOS,1.0D+00,
     &           SCR1,NMOS,SCR2,NMOS,0.0D+00,SCR3,NMOS)
      CALL XGEMM('N','N',NMOS,NMOS,NMOS,1.0D+00,
     &           SCR2,NMOS,SCR3,NMOS,0.0D+00,SCR1,NMOS)
C
C     Replace original MO coefficients by orthonormal set, C*X**-1/2.
C
ccc      CALL GETBLK(C     ,SCR2,NBFIRR(I),NBAS,IREPS(I))
      CALL ICOPY(NBFIRR(IRREP)*NBFIRR(IRREP)*IINTFP,C,1,SCR2,1)
      CALL  XGEMM('N','N',NBFIRR(IRREP),NMOS,NMOS,1.0D+00,
     &            SCR2,NBFIRR(IRREP),SCR1,NMOS,0.0D+00,
     &            SCR3,NBFIRR(IRREP))
C
ccc      CALL PUTBLK(C,SCR3,NBFIRR(I),NBAS,IREPS(I))
c     CALL ICOPY(SCR3,C,NBFIRR(IRREP)*NBFIRR(IRREP)*IINTFP)
      CALL ICOPY(NBFIRR(IRREP)*NBFIRR(IRREP)*IINTFP,SCR3,1,C,1)
C
C     Compute X now with new set. C is in SCR3.
C
cc      CALL GETBLK(SOVRLP,SCR2,NBFIRR(IRREP),NBAS,IREPS(IRREP))
cc      CALL XGEMM('N','N',NBFIRR(I),NBFIRR(I),NBFIRR(I), 1.0D+00,
cc     &           SCR2,NBFIRR(I),SCR3,NBFIRR(I), 0.0D+00,
cc     &           SCR1,NBFIRR(I))
cc      CALL XGEMM('T','N',NBFIRR(I),NBFIRR(I),NBFIRR(I), 1.0D+00,
cc     &           SCR3,NBFIRR(I),SCR1,NBFIRR(I), 0.0D+00,
cc     &           SCR2,NBFIRR(I))
C
cc      IF(IFLAGS(1).GE.5)THEN
cc       CALL OUTPUT(SCR2,1,NBFIRR(I),1,NBFIRR(I),NBFIRR(I),NBFIRR(I),1)
cc      ENDIF
C                 
cc  100 CONTINUE
      RETURN
      END
