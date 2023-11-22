      SUBROUTINE R000(JMAX,NOINT,ISTEPT,ISTEPU,ISTEPV,NRTUV,NUCAB,NUCCD,
     *                THRESH,ONECEN,IPRINT)
C
C     TUH 84
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
      PARAMETER (ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0)
      LOGICAL ONECEN, NOINT
      COMMON /PRIVEC/ NPRM12(MXAOSQ), NPRM34(MXAOSQ),
     *                LSTP12(MXAOSQ), LSTP34(MXAOSQ),
     *                LSTO12(MXAOSQ), LSTO34(MXAOSQ),
     *                LASTP(MXAOSQ),  LASTO(MXAOSQ)
      COMMON /ODCCOM/ COORPX(MXAOSQ), COORQX(MXAOSQ), COORPY(MXAOSQ),
     *                COORQY(MXAOSQ), COORPZ(MXAOSQ), COORQZ(MXAOSQ),
     *                EXP12(MXAOSQ),  EXP34(MXAOSQ),
     *                FAC12(MXAOSQ),  FAC34(MXAOSQ)
C     GAMCOM: 3267 = 27*121, max J value = 20
C
      PARAMETER (MXQN=8)
      PARAMETER (MXAQN=MXQN*(MXQN+1)/2,MXAQNS=MXAQN*MXAQN)
      PARAMETER (MAXJ = 32)
      COMMON/GAMCOM/WVAL,FJW(0:4*(MXQN-1)+2),
     *              TABFJW(121*(4*(MXQN-1)+2+7)), JMAX0
C
#include "lwork.par"
      PARAMETER (MAXVEC = LWORK3 / 6)
      COMMON /CWORK3/ WK3LOW, DPQX(MAXVEC), DPQY(MAXVEC), DPQZ(MAXVEC),
     *                RJ000(3*MAXVEC)
      LOGICAL TKTIME
      LOGICAL         DC101, DC1H1, DC1E1, DC2H1, DC2E1,
     *                DC102, DC1H2, DC1E2, DC2H2, DC2E2,
     *                DPATH1, DPATH2
      COMMON /SUBDIR/ DC101, DC1H1, DC1E1, DC2H1, DC2E1,
     *                DC102, DC1H2, DC1E2, DC2H2, DC2E2,
     *                DPATH1, DPATH2, NTOTAL
      IF (.NOT.DPATH1) THEN
         NUC12  = NUCCD
         NUC34  = NUCAB
         IOFF12 = MXAOSQ
         IOFF34 = 0
         SIGN   = - ONE
      ELSE
         NUC12  = NUCAB
         NUC34  = NUCCD
         IOFF12 = 0
         IOFF34 = MXAOSQ
         SIGN   = ONE
      END IF
      IF (ISTEPT .GT. MAXVEC) THEN
         WRITE (LUPRI,'(2(A,I7),A)') ' ISTEPT = ', ISTEPT,
     *    ' GREATER THAN ',MAXVEC,' IN R000 - PROGRAM CANNOT PROCEED. '
         WRITE (LUPRI,'(/A)')
     *      ' POSSIBLE ACTIONS: DECREASE ".MAXPRI" IN "*READIN" INPUT',
     *      '                   OR INCREASE SPACE IN /CWORK3/.'
         CALL ERREX 
      END IF
      JMSTP = ISTEPT*(JMAX + 1)
      IF (JMSTP .GT. 3*MAXVEC) THEN
         WRITE (LUPRI,'(2(A,I7),A)') ' JMSTP = ',JMSTP,' GREATER THAN ',
     *     3*MAXVEC, ' IN R000 - PROGRAM CANNOT PROCEED. '
         WRITE (LUPRI,'(/A)')
     *      ' POSSIBLE ACTIONS: DECREASE ".MAXPRI" IN "*READIN" INPUT',
     *      '                   OR INCREASE SPACE IN /CWORK3/.'
         CALL ERREX
      END IF
      JMAX0 = JMAX
      NOINT = .TRUE.
      DO 110 I = 1, NUC12
         LASTP(I) = NUC34
  110 CONTINUE
C
C     *****************
C     ***** RJ000 *****
C     *****************
C
C     SPECIAL CASE: ONE-CENTER INTEGRALS
C
      IF (ONECEN) THEN
         IODS = 0
         IADR12 = IOFF12
         DO 200 IOD12 = 1, NUC12
            IADR12 = IADR12 + 1
            EXPP = EXP12(IADR12)
            FAC12I = FAC12(IADR12)
            IADR34 = IOFF34
            DO 210 IOD34 = 1, LASTP(IOD12)
               IADR34 = IADR34 + 1
               IODS = IODS + 1
               EXPQ = EXP12(IADR34)
               FAC34I = FAC12(IADR34)
               EXPPQ = EXPP + EXPQ
               FACTOR = FAC12I*FAC34I/SQRT(EXPPQ)
               IF (ABS(FACTOR) .GT. THRESH) THEN
                  NOINT = .FALSE.
                  ALPHA = EXPP*EXPQ/EXPPQ
                  DPQX(IODS) = ZERO
                  DPQY(IODS) = ZERO
                  DPQZ(IODS) = ZERO
                  TALPHA = ALPHA + ALPHA
                  IADR = IODS
                  FJ0INV = ONE
                  DO 220 I = 0, JMAX
                     RJ000(IADR) = FACTOR/FJ0INV
                     FACTOR = - TALPHA*FACTOR
                     FJ0INV = FJ0INV + TWO
                     IADR = IADR + ISTEPT
  220             CONTINUE
               ELSE
                  IADR = IODS
                  DO 225 I = 0 ,JMAX
                     RJ000(IADR) = ZERO
                     IADR = IADR + ISTEPT
  225             CONTINUE
               END IF
  210       CONTINUE
  200    CONTINUE
C
C     GENERAL CASE: MULTICENTER INTEGRALS
C
      ELSE
         IODS = 0
         IADR12 = IOFF12
         DO 300 IOD12 = 1, NUC12
            IADR12 = IADR12 + 1
            EXPP = EXP12(IADR12)
            FAC12I = FAC12(IADR12)
            PX = COORPX(IADR12)
            PY = COORPY(IADR12)
            PZ = COORPZ(IADR12)
            IADR34 = IOFF34
            DO 310 IOD34 = 1, LASTP(IOD12)
               IADR34 = IADR34 + 1
               IODS = IODS + 1
               EXPQ = EXP12(IADR34)
               FAC34I = FAC12(IADR34)
               EXPPQ = EXPP + EXPQ
               FACTOR = FAC12I*FAC34I/SQRT(EXPPQ)
               IF (ABS(FACTOR) .GT. THRESH) THEN
                  NOINT = .FALSE.
                  ALPHA = EXPP*EXPQ/EXPPQ
                  PQX = PX - COORPX(IADR34)
                  PQY = PY - COORPY(IADR34)
                  PQZ = PZ - COORPZ(IADR34)
                  DPQX(IODS) = SIGN*PQX
                  DPQY(IODS) = SIGN*PQY
                  DPQZ(IODS) = SIGN*PQZ
                  WVAL = ALPHA*(PQX*PQX + PQY*PQY + PQZ*PQZ)
                  CALL GAMFUN
                  TALPHA = ALPHA + ALPHA
                  IADR = IODS
                  DO 320 I = 0, JMAX
                     RJ000(IADR) = FACTOR*FJW(I)
                     FACTOR = - TALPHA*FACTOR
                     IADR = IADR + ISTEPT
  320             CONTINUE
               ELSE
                  IADR = IODS
                  DO 325 I = 0 ,JMAX
                     RJ000(IADR) = ZERO
                     IADR = IADR + ISTEPT
  325             CONTINUE
               END IF
  310       CONTINUE
  300    CONTINUE
      END IF
      IPRINT = 0
      IF (IPRINT .LT. 10) RETURN
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      WRITE (LUPRI, 1000)
      WRITE (LUPRI, 1010) JMAX
      WRITE (LUPRI, 1020) NUC12, NUC34
      WRITE (LUPRI, 1030) THRESH
      WRITE (LUPRI, 1040) ONECEN
      WRITE (LUPRI, 1045) NOINT
      WRITE (LUPRI, 1050) ISTEPT, ISTEPU, ISTEPV
      WRITE (LUPRI, 1060) NRTUV
 1000 FORMAT (//,' ********** SUBROUTINE R000 **********')
 1010 FORMAT (//,'  JMAX    ',I7)
 1020 FORMAT (   '  NUC     ',2I7)
 1030 FORMAT (   '  THRESH  ',1P,E12.4)
 1040 FORMAT (   '  ONECEN  ',L7)
 1045 FORMAT (   '  NOINT   ',L7)
 1050 FORMAT (   '  ISTEP   ',3I7)
 1060 FORMAT (   '  NRTUV   ',I7)
      IF (IPRINT .LT. 20) RETURN
      WRITE (LUPRI, 1100)
      ISTART = 0
      DO 2000 J = 0, JMAX
         WRITE (LUPRI, 1110) J
         IADR = ISTART
         DO 2100 I = 1, NUC12
            WRITE (LUPRI, 1120) I
            WRITE (LUPRI, 1130) (RJ000(IADR + K), K = 1, NUC34)
            IADR = IADR + NUC34
 2100    CONTINUE
         ISTART = ISTART + ISTEPT
 2000 CONTINUE
 1100 FORMAT(//,' ***** RJ000 - INTEGRALS *****')
 1110 FORMAT( /,'  JMAX  ',I7)
 1120 FORMAT( /,'  NUC12 ',I7,/)
 1130 FORMAT(1P,6E12.4)
      IPRINT =0
      RETURN
      END
