      SUBROUTINE HERDIS(ISTEPT,JMAX,NRTUV,IPRINT,WORK1,LWORK1)
C
C     TUH 87
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "ninahh.par"
      COMMON /CINAHH/ INAHH(NINAHH)
      COMMON /INTADR/ IWKAO, IWKSO, IWKHHS, IWK1HH, IWK1HC, IWKLST
      DIMENSION WORK1(LWORK1)
C
      IF ((JMAX + 1)**3 .GT. NINAHH) THEN
         WRITE (LUPRI,'(/1X,A/A,I5,A)')
     *         ' Insufficent memory for INAHH in HERDIS.',
     *         ' Increase NINAHH to ', (JMAX + 1)**3,'.'
         STOP 
      END IF
C
      IOFF = 0
      ISTVMU = JMAX*(JMAX + 1)
      DO 100 J = 0, JMAX
         IADRJ = J
         DO 200 IT = J, 0, -1
            IADR = IADRJ
            JMINT = J - IT
            DO 300 IU = JMINT, 0, -1
               IV = JMINT - IU
               ISTR1 = IWKHHS + NRTUV + IOFF
               ISTR2 = IWKHHS + ISTEPT*IADR
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
*VOCL LOOP,NOVREC
CDIR$ IVDEP
               DO 400 I = 1, ISTEPT
                  WORK1(ISTR1 + I) = WORK1(ISTR2 + I)
  400          CONTINUE
               INAHH(IADR + 1) = IWK1HH + IOFF
               IOFF = IOFF + ISTEPT
               IADR = IADR + ISTVMU
  300       CONTINUE
            IADRJ = IADRJ + JMAX
  200    CONTINUE
  100 CONTINUE
C
      IOFF = IWKHHS + NRTUV
      NR = (JMAX + 1)*(JMAX + 2)*(JMAX + 3)/6
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
*VOCL LOOP,NOVREC
CDIR$ IVDEP
      DO 500 I = 1, ISTEPT*NR
         WORK1(IWK1HH + I) = WORK1(IOFF + I)
  500 CONTINUE
      IF (IPRINT .GE. 25) THEN
         CALL HEADER('OUTPUT FROM HERDIS',-1)
         WRITE (LUPRI,'(/1X,A)') ' Hermitian integrals '
         IOFF = IWK1HH
         DO 600 I = 1, NR
            WRITE (LUPRI,'(A,I5/)') ' NR ',I
            WRITE (LUPRI,'(6F12.8)') (WORK1(IOFF+J),J=1,ISTEPT)
            WRITE (LUPRI,'()')
            IOFF = IOFF + ISTEPT
  600    CONTINUE
         WRITE (LUPRI,'(/1X,A/)')
     *      ' Indirect addresses for Hermitian integrals.'
         IOFF = 1
         DO 700 I = 0, JMAX
            DO 710 J = 0, JMAX
               WRITE (LUPRI,'(1X,10I7)')
     *          (INAHH(IOFF + K), K = 0, JMAX)
               IOFF = IOFF + JMAX + 1
  710       CONTINUE
            WRITE (LUPRI,'()')
  700    CONTINUE
      END IF
      RETURN
      END
