      SUBROUTINE PCCD_FINISH(ICYCLE,SIDE)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*6 STAR
      CHARACTER*1 SIDE
      LOGICAL PCCD,CCD,LCCD
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /ENERGY/ ECORR(500,2)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /CALC/PCCD,CCD,LCCD
      
      DATA TWO /2.0D0/

      CALL GETREC(20,'JOBARC','SCFENEG ',IINTFP,ESCF)
      IF (PCCD) THEN
      WRITE(6,99)
99    FORMAT(T3,'     Summary of iterative solution of PCCD equations ')
      ELSEIF (CCD) THEN
      WRITE(6,"(2a)") "     Summary of iterative solution of CCD", 
     &                " equations "
      ELSEIF (LCCD) THEN
      WRITE(6,"(2a)") "     Summary of iterative solution of LCCD", 
     &                " equations "
      ENDIF 
      WRITE(6,100)
100   FORMAT(6X,59('-'))
      WRITE(6,101)
101   FORMAT(T24,'Correlation',T46,'Total',/,T8,'Iteration',T26,
     &       'Energy',T46,'Energy')
      WRITE(6,100)

      DO 20 I=1,ICYCLE
       STAR = ' DIIS '
       IF (pCCD) THEN
           WRITE(6,1000)I-1,TWO*ECORR(I,1),TWO*ECORR(I,1)+ESCF,STAR
       ELSE
           WRITE(6,1000)I-1,TWO*ECORR(I,1),TWO*ECORR(I,1)+ESCF,STAR
       ENDIF 
1000   FORMAT(T10,I4,T19,F18.12,T39,F19.12,T59,A6)
20    CONTINUE
      WRITE(6,100)

      IF (SIDE .EQ. "R") THEN
      IF (PCCD) THEN
      WRITE(6,202)
      ELSE IF (CCD) THEN
      WRITE(6,203)
      ELSE IF (LCCD) THEN
      WRITE(6,204)
      ENDIF 
      ELSEIF (SIDE .EQ. "L") THEN
      IF (PCCD) THEN
      WRITE(6,205)
      ELSE IF (CCD) THEN
      WRITE(6,206)
      ELSE IF (LCCD) THEN
      WRITE(6,207)
      ENDIF 
      ENDIF 

      CALL PUTREC(1,'JOBARC','TOTENERG',IINTFP,ECORR(ICYCLE,1)+ESCF)
      IF (SIDE .EQ. "R") THEN
202   FORMAT(T7,'A miracle has come to pass. ',
     &          'The right PCCD iterations have converged.')
203   FORMAT(T7,'A miracle has come to pass. ',
     &          'The right CCD iterations have converged.')
204   FORMAT(T7,'A miracle has come to pass. ',
     &          'The right LCCD iterations have converged.')
      IF (PCCD) THEN
      WRITE(6,1010) "pCCD ",TWO*ECORR(ICYCLE,1)+ESCF
      ELSE IF (CCD) THEN
      WRITE(6,1010) "CCD ", TWO*ECORR(ICYCLE,1)+ESCF
      ELSE IF (LCCD) THEN
      WRITE(6,1010) "LCCD ",TWO*ECORR(ICYCLE,1)+ESCF
      ENDIF 
      ELSEIF (SIDE .EQ. "L") THEN
205   FORMAT(T7,'A miracle has come to pass. ',
     &          'The left PCCD iterations have converged.')
206   FORMAT(T7,'A miracle has come to pass. ',
     &          'The left CCD iterations have converged.')
207   FORMAT(T7,'A miracle has come to pass. ',
     &          'The left LCCD iterations have converged.')
      IF (PCCD) THEN
      WRITE(6,1020) "pCCD ",TWO*ECORR(ICYCLE,1)+ESCF
      ELSE IF (CCD) THEN
      WRITE(6,1020) "CCD ", TWO*ECORR(ICYCLE,1)+ESCF
      ELSE IF (LCCD) THEN
      WRITE(6,1020) "LCCD ",TWO*ECORR(ICYCLE,1)+ESCF
      ENDIF 
      ENDIF 
  
1010  FORMAT(/,T14,A13,' energy is ',F20.12,' a.u. ',/)
1020  FORMAT(/,T14,A13,' pseudo energy is ',F20.12,' a.u. ',/)

      RETURN
      END
