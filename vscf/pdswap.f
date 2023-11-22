      SUBROUTINE PDSWAP(EVEC,IANGBF,SCR,NBAS,iblk)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXANG = 4)
C
      LOGICAL I2MANYP,I2MANYD
C
      DIMENSION EVEC(NBAS,NBAS),SCR(NBAS,NBAS),IANGBF(NBAS)
      DIMENSION ISIZ(MAXANG)
C
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA IONE /1/
      DATA ISIZ /3,6,10,15/
C
C  If this is a generally contracted basis set, then we need to reorder
C  rows of the eigenvector matrix before doing the symmetry analysis.
C
      I2MANYP=.FALSE.
      I2MANYD=.FALSE.
      CALL ZERO(SCR,NBAS*NBAS)
      ICNT=0
      IJUMP=0
      DO 100 I=1,NBAS
        IF(I.LE.IJUMP) GOTO 100
        IF(IANGBF(I).GE.1) THEN
          IQ=ISIZ(IANGBF(I))
 3000     ICNT=ICNT+1
          IF(IANGBF(I+IQ*ICNT).EQ.1.AND.(I+IQ*ICNT).LE.NBAS) GOTO 3000
          IF(ICNT.GT.7) I2MANYP=.TRUE.
          DO 110 J=I,I+IQ*ICNT-1
            INEW=MOD(J-I,ICNT)*IQ+((J-I)/ICNT)+I
            CALL SCOPY(NBAS,EVEC(J,1),NBAS,SCR(INEW,1),NBAS)
  110     CONTINUE
          IJUMP=I+IQ*ICNT-1
          ICNT=0
        ELSE
          CALL SCOPY(NBAS,EVEC(I,1),NBAS,SCR(I,1),NBAS)
        ENDIF
  100 CONTINUE
c YAU : old
c     CALL ICOPY(NBAS*NBAS*IINTFP,SCR,1,EVEC,1)
c YAU : new
      CALL DCOPY(NBAS*NBAS,SCR,1,EVEC,1)
c YAU : end
C
      IF((I2MANYP.OR.I2MANYD).and.iblk.eq.1) THEN
       WRITE(LUOUT,5000)
 5000  FORMAT(/,T3,'@PDSWAP-W,',5X,20('*'),'   WARNING   ',
     &             20('*'),/,/,
     &        T20,'There may be more than 7 p-functions or more than',/,
     &        T20,'7 d-functions per center in this calculation.',/,
     &        T20,'This causes a new shell to be used to include',/,
     &        T20,'the additional functions.  This is outside the',/,
     &        T20,'primary assumption of PDSWAP.  The results MAY',/,
     &        T20,'be incorrect.',/,/,T18,53('*'),/)
      ENDIF
      RETURN
      END
