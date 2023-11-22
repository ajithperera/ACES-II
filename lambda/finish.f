C
C RUN-DOWN ROUTINE FOR CC CALCULATIONS.
C
      SUBROUTINE FINISH(ICYCLE,pCCD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*6 STAR
      LOGICAL pCCD
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
cAP - 500 is the maximum number of CC iterations (not basis functions)
      COMMON /ENERGY/ ECORR(500,2),IXTRLE(500)
      COMMON /FLAGS/ IFLAGS(100)
C
CJDW 10/1/96. Stop this routine from writing out the summary of
C             lambda iterations, unless the print flag is set high.
C             If the summary is to be printed, emphasize the numbers
C             are pseudo-energies and are not the real energy. Now
C             no one should ever read the lambda energy as the CC
C             energy again. If anyone does, he/she will be asked to
C             rewrite aces2 in Fortran 66 and make it run on SGI, DEC,
C             and PC.
C
CSSS      IF(IFLAGS(1).GE.10)THEN
C
      DATA HALF /0.50D0/
      CALL GETREC(20,'JOBARC','SCFENEG ',IINTFP,ESCF)
      WRITE(6,99)
99    FORMAT(T3,
     & '     Summary of iterative solution of Lambda equations ')
      WRITE(6,100)
100   FORMAT(6X,59('-'))
      WRITE(6,101)
101   FORMAT(T26,'Lambda',T46,'Total',/,T8,'Iteration',T23,
     &       'Pseudo-energy',T43,'Pseudo-energy')
      WRITE(6,100)
      DO 20 I=1,ICYCLE
       STAR='JACOBI'
       IF(IXTRLE(I).EQ.1)STAR=' RLE  '
       IF(IFLAGS(21).EQ.1)STAR=' DIIS '
       IF (pCCD) THEN 
          WRITE(6,1000)I-1,ECORR(I,1),
     &                     ECORR(I,1)+ESCF,STAR
       ELSE
          WRITE(6,1000)I-1,ECORR(I,1),ECORR(I,1)+ESCF,STAR
       ENDIF 
1000   FORMAT(T10,I4,T19,F18.12,T39,F19.12,T59,A6)
20    CONTINUE
      WRITE(6,100)
      WRITE(6,201)
  201 FORMAT(/,6X,59('-'),/,6X,
     & ' ** Warning : The Pseudo-energy is NOT the real energy. **',
     & /,6X,
     & '          The real energy is printed out by xvcc          ',
     & /,6X,59('-'))
      WRITE(6,202)
202   FORMAT(T7,'A miracle has come to pass. ',
     &          'The Lambda iterations have converged.')
C
CSSS      ENDIF
C
      WRITE(6,1010)
 1010 FORMAT(/,77('-'),/,T30,'Exiting xlambda',/,77('-'))

C For pCCD, restore the CCD for any module beyond lambda, 10/2021

      IF (pCCD) THEN 
         IFLAGS(h_IFLAGS_calc) = 8
         CALL PUTREC(20,"JOBARC","PCCD_RUN",1,1)
         CALL PUTREC(20,"JOBARC","IFLAGS  ",100,IFLAGS)
      ENDIF 
      call aces_fin
      STOP
      RETURN
      END
