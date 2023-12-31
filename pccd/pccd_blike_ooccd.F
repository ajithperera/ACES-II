      SUBROUTINE Pccd_blike_ooccd(ICORE,MAXCOR,IUHF)
C
C CONTROL MODULE FOR BRUECKNER CALCULATIONS
C
      IMPLICIT INTEGER (A-H,O-Z)
      DOUBLE PRECISION TOL,X,FNDLRGAB
      DIMENSION ICORE(MAXCOR)
      LOGICAL BRKCNV, BRUKMOS_EXIST
      LOGICAL pCCD
      COMMON  /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /FLAGS/ IFLAGS(100)
#include "info.com"
C
C GET TOLERANCE FROM IFLAG ARRAY
C
      IONE=1
      TOL=10.0**(-IFLAGS(76))
C
C CHECK CURRENT T1 AMPLITUDES TO SEE IF CONVERGENCE HAS BEEN
C  ACHIEVED
C
      I000=1
      I010=I000+IINTFP*NT(1)
      IF(IUHF.NE.0)I020=I010+IINTFP*NT(2)
      LEN=NT(1)+IUHF*NT(2)
      CALL GETLST(ICORE(I000),1,1,1,1,90)
      IF(IUHF.NE.0)CALL GETLST(ICORE(I010),1,1,1,2,90)
      X=FNDLRGAB(ICORE(I000),LEN)
      WRITE(6,90)
90    FORMAT(T2,'Determination of optimized OO-CCD  orbitals:') 
      WRITE(6,100)X
100   FORMAT(T2,'Largest orbital rotation gradeint is ',F20.10,'.')
      IF(X.LT.TOL)THEN
       CALL PUTREC(20,'JOBARC','BRUKTEST',IONE,1)
       WRITE(6,200)
200    FORMAT(T2,'OO-LCCD/OO-OCCD/MBPT(2) calculation converged.')
       BRKCNV = .TRUE.
       NBAS = NOCCO(1) + NVRTO(1)

      ELSE

       WRITE(6,300)
300    FORMAT(T3,"Orbital Rotation gradients are used to determine new",
     +           " orbitals.")
       CALL PUTREC(20,'JOBARC','BRUKTEST',IONE,0)
       CALL PCCD_ROTT1(ICORE,MAXCOR,IUHF)
       IF(IFLAGS(11).GT.0)IFLAGS(11)=1
       IFLAGS(16)=0
       IFLAGS(38)=1
       IFLAGS(34)=0
       BRKCNV = .FALSE.
      ENDIF
      CALL PUTREC(20,'JOBARC','IFLAGS  ',100,IFLAGS)
      RETURN
      END
