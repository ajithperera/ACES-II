      SUBROUTINE DENW101(EVAL,T1,IUHF,NBAS)
C
C DIVIDE SINGLE AMPLITUDES BY ENERGY DENOMINATORS AND WRITE
C TO WORKING AMPLITUDE LISTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT,A,B
      DIMENSION T1(1),EVAL(1)
      CHARACTER*8 LABEL(2)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPX(255,2),DIRPRD(8,8)
      DATA LABEL /'SCFEVALA','SCFEVALB'/
      DATA ZILCH /0.0/
C
C READ IN EIGENVALUES
C
      DO 10 ISPIN=1,1+IUHF
       CALL GETREC(20,'JOBARC',LABEL(ISPIN),IINTFP*NBAS,EVAL)
       CALL GETLST(T1,1,1,1,2+ISPIN,94)
       IOFF=0
       ITHRU=1
       DO 20 IRREP=1,NIRREP
        NOCC=POP(IRREP,ISPIN)
        DO 30 I=1,NOCC
         DO 31 J=1,NOCC
          Z=EVAL(IOFF+I)-EVAL(IOFF+J)
          IF(Z.NE.ZILCH)T1(ITHRU)=T1(ITHRU)/Z
          ITHRU=ITHRU+1
31       CONTINUE
30      CONTINUE
        IOFF=IOFF+NOCC
20     CONTINUE
       CALL PUTLST(T1,1,1,1,ISPIN,94)
10    CONTINUE
C
      RETURN
      END
