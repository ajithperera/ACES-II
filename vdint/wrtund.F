      SUBROUTINE WRTUND(NBAST,KSMAT,KHMAT,KTMAT,IPRINT,WORK1,LWORK1)
C
C     T. Helgaker Sep. 1986
C
C     This subroutine writes undifferentiated overlap, one-electron
C     Hamiltonian, and kinetic energy integrals on sequential file.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      COMMON /HRUNIT/ LUONEL, LUINTA, LUSUPM, LUDASP, LUSOL
C
      DIMENSION WORK1(LWORK1)
#include "lwork.par"
      COMMON /CWORK2/ WK2LOW, WORK2(LWORK2), WK2HGH
      COMMON /CCOM/ THRS, NHTYP, IBX
C
      DIMENSION   BUF(600), IBUF(600)
      EQUIVALENCE(WORK2(1),BUF(1)), (WORK2(601),IBUF(1))
      CHARACTER*8 LAB123(3), TITLES, TITLEH, TITLET
      DATA LAB123/'********','        ','        '/
      DATA TITLES/'OVERLAP '/, TITLEH/'ONEHAMIL'/, TITLET/'KINETINT'/
      ITRI(I,J) = I*(I-1)/2 + J
      INDEX1(IJ) = INT(SQRT(FLOAT(2*IJ) + 0.25) + 0.4999)
      INDEX2(IJ) = IJ - INDEX1(IJ)*(INDEX1(IJ)-1)/2
C
      IF (IPRINT .GE. 5) CALL HEADER ('OUTPUT FROM WRTUND',-1)
      DO 100 ITYPE = 1, 3
         IF (ITYPE .EQ. 1) THEN
            IF (IPRINT .GE. 5) WRITE (LUPRI,'(/A)') ' Overlap matrix'
            WRITE (LUONEL) LAB123, TITLES
            IOFF = KSMAT - 1
         ELSE IF (ITYPE .EQ. 2) THEN
            IF (IPRINT .GE. 5)
     *         WRITE (LUPRI,'(/A)')' One-el. Ham. matrix'
            WRITE (LUONEL) LAB123, TITLEH
            IOFF = KHMAT - 1
         ELSE IF (ITYPE .EQ. 3) THEN
            IF (IPRINT .GE. 5) WRITE (LUPRI,'(/A)')' Kin. energy matrix'
            WRITE (LUONEL) LAB123, TITLET
            IOFF = KTMAT - 1
         END IF
         IAB = 0
         ICOUNT = 0
         NBUF = 0
         DO 200 IA = 1, NBAST
            DO 300 IB = 1, IA
               IAB = IAB + 1
               AINT = WORK1(IOFF + IAB)
               IF (ABS(AINT) .GT. THRS) THEN
                  ICOUNT = ICOUNT + 1
                  BUF(ICOUNT) = AINT
                  IBUF(ICOUNT) = ITRI(IA,IB)
                  IF (ICOUNT .EQ. 600) THEN
                     WRITE (LUONEL) BUF, IBUF, ICOUNT
                     NBUF = NBUF + 1
                     IF (IPRINT .GE. 5) THEN
                        WRITE (LUPRI, '(//,A,2I5,/)')
     *                     ' ICOUNT, NBUF ', ICOUNT, NBUF
                        DO 400 I = 1, ICOUNT
                           LABEL = IBUF(I)
                           IORBA = INDEX1(LABEL)
                           IORBB = INDEX2(LABEL)
                           WRITE (LUPRI,'(1X,3I5,5X,E24.12)')
     *                           LABEL, IORBA, IORBB, BUF(I)
  400                   CONTINUE
                     END IF
                     ICOUNT = 0
                  END IF
               END IF
  300       CONTINUE
  200    CONTINUE
         NCOUNT = 600*NBUF + ICOUNT
         IF (ICOUNT .GT. 0) THEN
            WRITE (LUONEL) BUF, IBUF, ICOUNT
            NBUF = NBUF + 1
         END IF
         IF (IPRINT .GE. 5) THEN
             WRITE (LUPRI, '(//,A,2I5,/)')
     *             ' ICOUNT, NBUF ', ICOUNT, NBUF
             DO 500 I = 1, ICOUNT
                LABEL = IBUF(I)
                IORBA = INDEX1(LABEL)
                IORBB = INDEX2(LABEL)
                WRITE (LUPRI,'(1X,3I5,5X,E24.12)')
     *                LABEL, IORBA, IORBB, BUF(I)
  500        CONTINUE
         END IF
         WRITE (LUONEL) BUF, IBUF, -1
         PERCNT = FLOAT(100*NCOUNT) / FLOAT(IAB)
         IF (ITYPE .EQ. 1) THEN
            WRITE (LUPRI,'(//,I6,A,I4,A)') NCOUNT,
     *        ' atomic overlap integrals written in',NBUF,' buffers.'
         ELSE IF (ITYPE .EQ. 2) THEN
            WRITE (LUPRI,'(//,I6,A,I4,A)') NCOUNT,
     *        ' one-el. Hamil. integrals written in',NBUF,' buffers.'
         ELSE
            WRITE (LUPRI,'(//,I6,A,I4,A)') NCOUNT,
     *        ' kinetic energy integrals written in',NBUF,' buffers.'
         END IF
         WRITE(LUPRI, '(A,F6.2)') ' Percentage non-zero integrals: ',
     *                              PERCNT
  100 CONTINUE
      RETURN
      END
