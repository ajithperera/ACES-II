      SUBROUTINE BRUECKIT(ICORE,MAXCOR,IUHF,BRKCNV,pCCD)
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
90    FORMAT(T3,'Determination of Brueckner orbitals:') 
      WRITE(6,100)X
100   FORMAT(T3,'Largest T1 amplitude is ',F20.10,
     &          '.')
      IF(X.LT.TOL)THEN
       CALL PUTREC(20,'JOBARC','BRUKTEST',IONE,1)
       WRITE(6,200)
200    FORMAT(T3,'Brueckner CC/MBPT calculation converged.')
       BRKCNV = .TRUE.
C      
       NBAS = NOCCO(1) + NVRTO(1)
C This is incorrect for symmetry on. Can be fixed but the NEWMOS
C files has the same info. 

       CALL WRITE_BRUECKMOS(ICORE(I000),MAXCOR,IUHF,NBAS)

      ELSE
       WRITE(6,300)
300    FORMAT(T3,'T1 amplitudes are used to determine new orbitals.')
       CALL PUTREC(20,'JOBARC','BRUKTEST',IONE,0)
       CALL ROTT1(ICORE,MAXCOR,IUHF)
       IF(IFLAGS(11).GT.0)IFLAGS(11)=1
       IFLAGS(16)=0
       IFLAGS(38)=1
       IFLAGS(34)=0
C IFLAGS(77) is QRHF_GEN. I am not quite sure now why I commented 
C this. This is done during S-T gap of biradical studies. Ajith Perera
C 12/2013. 
CSSS       IFLAGS(77)=0
       BRKCNV = .FALSE.
C IF pCCD, change the calclevel back to pCCD from CCD.
      IF (pCCD) IFLAGS(2) = 52
      ENDIF
      CALL PUTREC(20,'JOBARC','IFLAGS  ',100,IFLAGS)
      RETURN
      END
