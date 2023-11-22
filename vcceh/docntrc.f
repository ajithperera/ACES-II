C
      SUBROUTINE DOCNTRC(Q, W, T, DISSYQ, NUMSYQ, DISSYT, NUMSYT,
     &                   DISSYW, NUMSYW, IRREPR, LISTW, LISTW0, LISTT,
     &                   CALMOD, C0, COEFF, ISWTCH)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*3 CALMOD
      INTEGER DISSYT,DISSYW,DISSYQ
      DIMENSION W(DISSYW,NUMSYW),T(DISSYT,NUMSYT),Q(DISSYQ,NUMSYQ)
      DATA TWO /2.0/
      DATA ONEM /-1.0/
C
C Pick up relevant integrals.
C     
      IF(ISWTCH .EQ. 0) THEN
C     
         CALL GETLST(W, 1, NUMSYW, 2, IRREPR, LISTW)
         IF(CALMOD .EQ. 'TxW')THEN
            CALL XGEMM('N', 'N', DISSYT, NUMSYW, NUMSYT, COEFF, T, 
     &                  DISSYT, W, DISSYW, C0, Q, DISSYQ)
         ELSE IF (CALMOD .EQ. 'WxT') THEN
            CALL XGEMM('N', 'N', DISSYW, NUMSYT, NUMSYW, COEFF, W,
     &                  DISSYW, T, DISSYT, C0, Q, DISSYQ)
         ENDIF
C     
      ELSE
C     
C Spin adapted code.
C     
         CALL GETLST(W, 1, NUMSYW, 2, IRREPR, LISTW)
         CALL GETLST(Q, 1, NUMSYW, 2, IRREPR, LISTW0)
C
         CALL SSCAL(NUMSYW*NUMSYW, TWO, W, 1)
         CALL SAXPY(NUMSYW*NUMSYW, ONEM,  Q, 1, W, 1)
         IF(CALMOD .EQ. 'TxW') THEN
            CALL XGEMM('N', 'N', DISSYT, NUMSYW, NUMSYT, COEFF, T,
     &                  DISSYT, W, DISSYW, C0, Q, DISSYQ)
         ELSEIF(CALMOD .EQ. 'WxT') THEN
            CALL XGEMM('N', 'N', DISSYW, NUMSYT, NUMSYW, COEFF, W,
     &                  DISSYW, T, DISSYT, C0, Q, DISSYQ)
         ENDIF
C     
      ENDIF
C
      RETURN
      END
