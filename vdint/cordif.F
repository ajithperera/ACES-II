      SUBROUTINE CORDIF(NORBA,NORBB,THRESH,D0X,D0Y,D0Z,
     *                  IPRINT,NUCOA,NUCOB,JSTRA,JSTRB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "baslims.par"
      LOGICAL DIFX, DIFY, DIFZ, D0X, D0Y, D0Z
      DIMENSION NUCOA(MXAOVC), NUCOB(MXAOVC),
     *          JSTRA(MXAOVC), JSTRB(MXAOVC)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      COMMON /PRIMIT/ PRIEXP(MXPRIM), PRICCF(MXPRIM,MXCONT),
     *                PRICRX(MXPRIM), PRICRY(MXPRIM), PRICRZ(MXPRIM)
C
      DIFX = .FALSE.
      DIFY = .FALSE.
      DIFZ = .FALSE.
C
C     A - A
C
      DO 100 I = 1, NORBA
      DO 100 J = 1, NUCOA(I)
         IJ  = JSTRA(I) + J
         CRX = PRICRX(IJ)
         CRY = PRICRY(IJ)
         CRZ = PRICRZ(IJ)
         DO 200 K = 1, I
         DO 200 L = 1, NUCOA(K)
            KL   = JSTRA(K) + L
            DIFX = DIFX .OR. ABS(PRICRX(KL)-CRX) .GT. THRESH
            DIFY = DIFY .OR. ABS(PRICRY(KL)-CRY) .GT. THRESH
            DIFZ = DIFZ .OR. ABS(PRICRZ(KL)-CRZ) .GT. THRESH
  200    CONTINUE
  100 CONTINUE
C
C     B - B
C
      IF (.NOT.(DIFX .AND. DIFY .AND. DIFZ)) THEN
         DO 300 I = 1, NORBB
         DO 300 J = 1, NUCOB(I)
            IJ  = JSTRB(I) + J
            CRX = PRICRX(IJ)
            CRY = PRICRY(IJ)
            CRZ = PRICRZ(IJ)
            DO 400 K = 1, I
            DO 400 L = 1, NUCOB(K)
               KL   = JSTRB(K) + L
               DIFX = DIFX .OR. ABS(PRICRX(KL)-CRX) .GT. THRESH
               DIFY = DIFY .OR. ABS(PRICRY(KL)-CRY) .GT. THRESH
               DIFZ = DIFZ .OR. ABS(PRICRZ(KL)-CRZ) .GT. THRESH
  400       CONTINUE
  300    CONTINUE
      END IF
C
C     A - B
C
      IF (.NOT.(DIFX .AND. DIFY .AND. DIFZ)) THEN
         DO 500 I = 1, NORBA
         DO 500 J = 1, NUCOA(I)
            IJ  = JSTRA(I) + J
            CRX = PRICRX(IJ)
            CRY = PRICRY(IJ)
            CRZ = PRICRZ(IJ)
            DO 600 K = 1, NORBB
            DO 600 L = 1, NUCOB(K)
               KL   = JSTRB(K) + L
               DIFX = DIFX .OR. ABS(PRICRX(KL)-CRX) .GT. THRESH
               DIFY = DIFY .OR. ABS(PRICRY(KL)-CRY) .GT. THRESH
               DIFZ = DIFZ .OR. ABS(PRICRZ(KL)-CRZ) .GT. THRESH
  600       CONTINUE
  500    CONTINUE
      END IF
      D0X = .NOT.DIFX
      D0Y = .NOT.DIFY
      D0Z = .NOT.DIFZ
C
      IF (IPRINT .LT. 05) RETURN
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      CALL HEADER('SUBROUTINE CORDIF',-1)
      WRITE (LUPRI, 1010) NORBA, NORBB
      WRITE (LUPRI, 1020) (NUCOA(I), I = 1, NORBA)
      WRITE (LUPRI, 1030) (NUCOB(I), I = 1, NORBB)
      WRITE (LUPRI, 1040) (JSTRA(I), I = 1, NORBA)
      WRITE (LUPRI, 1050) (JSTRB(I), I = 1, NORBB)
      WRITE (LUPRI, 1060) D0X, D0Y, D0Z
 1010 FORMAT(  '  NORB     ',2I7)
 1020 FORMAT(  '  NUCOA:   ',15I7)
 1030 FORMAT(  '  NUCOB:   ',15I7)
 1040 FORMAT(  '  JSTRA:   ',15I7)
 1050 FORMAT(  '  JSTRB:   ',15I7)
 1060 FORMAT(  '  D0X/Y/Z: ',3L5)
      RETURN
      END
