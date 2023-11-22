      SUBROUTINE COMPPR(PROP,DENS,PRPINT,NSIZ,NUCLEAR)
C
C THIS ROUTINE READS A LIST OF PROPERTY INTEGRALS AND THE DENSITY
C  MATRIX AND COMPUTES THE PARTICULAR PROPERTY.
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL NUCLEAR
      CHARACTER*32 CRAP
      DIMENSION DENS(NSIZ,NSIZ),PRPINT(NSIZ,NSIZ),BUF(600),IBUF(600)
      DIMENSION IXX(2)
      EQUIVALENCE (PRPNUC,IXX(1))
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      NNM1O2(IX)=(IX*(IX-1))/2
      IEXTI(IX)=1+(-1+INT(DSQRT(8.D0*IX+0.999D0)))/2
      IEXTJ(IX)=IX-NNM1O2(IEXTI(IX))
C
C FIRST GET NUCLEAR CONTRIBUTION BY BACKSPACING FILE.
C
      BACKSPACE(10)
      READ(10)CRAP,PRPNUC
      IF(.NOT.NUCLEAR) PRPNUC=0.D0
C
C READ IN THE PROPERTY INTEGRALS
C
      CALL ZERO(PRPINT,NSIZ*NSIZ)
1     READ(10)BUF,IBUF,NUT
CDIR$ IVDEP
*VOCL LOOP,NOVREC
      DO 10 I=1,NUT
       INDI=IEXTI(IBUF(I))
       INDJ=IEXTJ(IBUF(I))
       PRPINT(INDI,INDJ)=BUF(I)
       PRPINT(INDJ,INDI)=BUF(I)
10    CONTINUE
      IF(NUT.EQ.600)GOTO 1
C
C COMPUTE THE PROPERTY
C
      PROP=SDOT(NSIZ*NSIZ,DENS,1,PRPINT,1)+PRPNUC
      RETURN
      END
