      SUBROUTINE HERNAI(JMAX,PVAL,CPX,CPY,CPZ,FACINT,IADR,
     *                  ISTEPU,ISTEPV,NRTUV,IPRINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
C     This subroutine calculates the R integrals as defined by
C     McMurchie and Davidson in J. Comp. Phys. 26 (1978) 218.
C     The recursion formulas (4.6) - (4.8) are used.
C
C     Number of elements generated: JMAXP*(JMAXP+1)*(JMAXP+2)/6
C
C     The array RTUV(J) must be dimensioned as
C
C     J = (NUMNUC + 1)*((MAX(JMAX) + 1)**3)
C
C     where NUCNUM Is the total number of nuclei.
C
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
      PARAMETER (ONE = 1.D0, D2PI = 6.28318 53071 79586 D00)
#include "lwork.par"
#include "mxcent.par"
      PARAMETER (MXQN=8)
      PARAMETER (MXAQN=MXQN*(MXQN+1)/2,MXAQNS=MXAQN*MXAQN)
      PARAMETER (MAXJ = 32)
      PARAMETER (LWKRLM = LWORK3 - 8020)
      COMMON /CWORK3/ WK3LOW, SHGTF, RTUV((MXCENT+1)*(2*MXQN+1)**3)
CSSS                  RLMCOF(LWKRLM), WK3HGH
C     GAMCOM: 3267 = 27*121, max J value = 20
CSSS      PARAMETER (MAXJ = 20)
CSSS      COMMON /GAMCOM/ WVAL, FJW(0:MAXJ), TABFJW(3267), JMAX0
C
      COMMON/GAMCOM/WVAL,FJW(0:4*(MXQN-1)+2),
     *              TABFJW(121*(4*(MXQN-1)+2+7)), JMAX0

      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
      IF (IPRINT .GT. 20) THEN
         WRITE (LUPRI,'(/A/)') ' <<<<< Output from Hernai >>>>> '
         WRITE (LUPRI,'(A,I10)')   ' JMAX   ', JMAX
         WRITE (LUPRI,'(A,I10)')   ' IADR   ', IADR
         WRITE (LUPRI,'(A,I10)')   ' ISTEPU ', ISTEPU
         WRITE (LUPRI,'(A,I10)')   ' ISTEPV ', ISTEPV
         WRITE (LUPRI,'(A,I10)')   ' NRTUV  ', NRTUV
         WRITE (LUPRI,'(A,F12.6)') ' PVAL   ', PVAL
         WRITE (LUPRI,'(A,F12.6)') ' FACINT ', FACINT
         WRITE (LUPRI,'(A,3F12.6)') ' CP ', CPX, CPY, CPZ
      END IF
C
C  Incomplete Gamma Function 
C
      WVAL = PVAL*(CPX*CPX + CPY*CPY + CPZ*CPZ)
      JMAX0 = JMAX
      CALL GAMFUN
      IF (IPRINT .GT. 20) THEN
       CALL HEADER ('FJW after GAMFUN',-1)
       WRITE (LUPRI,'(6F20.10)') (FJW(I), I = 1, JMAX)
      END IF

C
C   Special Case: JMAX = 0 
C
      IF (JMAX .EQ. 0) THEN
       RTUV(IADR) = FACINT*D2PI*FJW(0)/PVAL
       RETURN
      END IF
C
C  General Case: JMAX > 0 
C
      IF (IBTAND(JMAX,1) .EQ. 0) THEN
       ISTRTJ = IADR
       ISTEPJ = NRTUV
      ELSE
       ISTRTJ = IADR + NRTUV
       ISTEPJ = - NRTUV
      END IF
      D2PVAL = PVAL + PVAL
      FACTOR = FACINT*D2PI/PVAL
      DO 100 JVAL = 0, JMAX
       FJW(JVAL) = FACTOR*FJW(JVAL)
       FACTOR    = - D2PVAL*FACTOR
  100 CONTINUE
      IF (IPRINT .GT. 20) THEN
       CALL HEADER ('FJW after multiplication by FACTOR',-1)
       WRITE (LUPRI,'(6F20.10)') (FJW(I), I = 1, JMAX)
      END IF
C
C  JVAL = 0 
C
      RTUV(ISTRTJ) = FJW(JMAX)
C
C  JVAL = 1 
C
      ISTRTJ                =   ISTRTJ + ISTEPJ
      FJWMAX                =   FJW(JMAX)
      RTUV(ISTRTJ)          =   FJW(JMAX - 1)
      RTUV(ISTRTJ +      1) = - CPX*FJWMAX
      RTUV(ISTRTJ + ISTEPU) = - CPY*FJWMAX
      RTUV(ISTRTJ + ISTEPV) = - CPZ*FJWMAX
C
      IF (JMAX .EQ. 1) RETURN
C
C  JVAL > 1 
C
      ISTPTU = 1 - ISTEPU
      ISTEPJ =   - ISTEPJ
      DO 200 JVAL = 2, JMAX
         ISTRTJ = ISTRTJ + ISTEPJ
         ISTRTV = ISTRTJ
         ISTRTU = ISTRTV
         IPREV1 = ISTRTU - ISTEPJ - 1
         IPREV2 = IPREV1 - 1
C
C        RJ(0,0,0)
C
         RTUV(ISTRTU) = FJW(JMAX - JVAL)
C
C        RJ(1,0,0)
C
         RTUV(ISTRTU + 1) = - CPX*RTUV(IPREV1 + 1)
C
C        RJ(t,0,0) for  t > 1
C
         TMIN1 = ONE
         DO 300 IT = 2, JVAL
            RTUV(ISTRTU + IT) = - CPX*RTUV(IPREV1 + IT)
     *                          + TMIN1*RTUV(IPREV2 + IT)
            TMIN1 = TMIN1 + ONE
  300    CONTINUE
C
C        RJ(t,1,0)
C
         ISTRTU = ISTRTU + ISTEPU
         IPREV1 = IPREV1 + 1
         DO 400 IT = 0, JVAL - 1
            RTUV(ISTRTU + IT) = - CPY*RTUV(IPREV1 + IT)
  400    CONTINUE
C
C        RJ(t,u,0) for  u > 1
C
         UMIN1 = ONE
         DO 500 IU = 2,JVAL
            ISTRTU = ISTRTU + ISTEPU
            IPREV1 = IPREV1 + ISTEPU
            IPREV2 = IPREV1 - ISTEPU
            DO 510 IT = 0, JVAL - IU
               RTUV(ISTRTU + IT) = - CPY*RTUV(IPREV1 + IT)
     *                             + UMIN1*RTUV(IPREV2 + IT)
  510       CONTINUE
            UMIN1 = UMIN1 + ONE
  500    CONTINUE
C
C        RJ(t,u,1)
C
         ISTRTV = ISTRTV + ISTEPV
         ISTRTU = ISTRTV
         IPREV1 = ISTRTU - ISTEPJ - ISTEPV
         IUMAX = JVAL - 1
         DO 600 IU = 0, IUMAX
            DO 610 IT = 0, IUMAX - IU
               RTUV(ISTRTU + IT) = - CPZ*RTUV(IPREV1 + IT)
  610       CONTINUE
            ISTRTU = ISTRTU + ISTEPU
            IPREV1 = IPREV1 + ISTEPU
  600    CONTINUE
C
C        RJ(t,u,v) for v > 1
C
         VMIN1 = ONE
         DO 700 IV = 2,JVAL
            ISTRTV = ISTRTV + ISTEPV
            ISTRTU = ISTRTV
            IPREV1 = ISTRTU - ISTEPJ - ISTEPV
            IPREV2 = IPREV1 - ISTEPV
            IUMAX = JVAL - IV
            DO 710 IU = 0, IUMAX
               DO 720 IT = 0, IUMAX - IU
                  RTUV(ISTRTU + IT) = - CPZ*RTUV(IPREV1 + IT)
     *                                + VMIN1*RTUV(IPREV2 + IT)
  720          CONTINUE
               ISTRTU = ISTRTU + ISTEPU
               IPREV1 = IPREV1 + ISTEPU
               IPREV2 = IPREV2 + ISTEPU
  710       CONTINUE
            VMIN1 = VMIN1 + ONE
  700    CONTINUE
         ISTEPJ = - ISTEPJ
  200 CONTINUE
      RETURN
      END
