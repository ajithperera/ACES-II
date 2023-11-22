      SUBROUTINE EXCOEF (TSF,L,M,NUCA,NUCB,JSTA,JSTB,FAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     Parameters
C-----------------------------------------------------------------------
#include "baslims.par"
#include "mxatms.par"
      parameter (nh4=4*NHT-3, nh2=NHT+NHT+1)
C-----------------------------------------------------------------------
      COMMON /DAT/  ALPHA(MXTNPR),CONT(MXTNCC),CENT(3,MXTNSH),
     1              CORD(Mxatms,3),CHARGE(Mxatms),FMULT(8),TLA, TLC
      DIMENSION TSF(40000),C(nh2),TSFA(nh2),TSFAB(nh2),FAC(nh4)
C
CJDW 9/9/97. I've tried altering this routine since it was not working
C            yesterday on crunch when optimization was used. The compil-
C            ation used was "f77 -c -dalign -fns -fast -fsimple=2 -O5".
C            I changed the Fortran so it looks like f77. The main problem
C            for the optimizer was probably the fact that calculation of
C            TSF and the update of TSFAB were in the same loop. Now I've
C            split them.
C
C            Note the 40000 limit.
C
C     IA   TSFA(1) TSFA(2) TSFA(3) TSFA(4) TSFA(5) TSFA(6) TSFA(7)
C
C      1     1
C      2     B        1
C      3     B**2    2B       1
C      4     B**3    3B**2   3B       1
C      5     B**4    4B**3   6B**2   4B       1
C      6     B**5    5B**4  10B**3  10B**2   5B       1
C      7     B**6    6B**5  15B**4  20B**3  15B**2   6B       1
C
C
      IADR=1
      LM=L+M-1
      C(1)=1.0
      DO 100 JA=1,NUCA
      A=ALPHA(JA+JSTA)
      DO  90 JB=1,NUCB
      B=-ALPHA(JB+JSTB)
      D=1./(A-B)
C
      DO 10 I=2,LM
      C(I)=C(I-1)*D
   10 CONTINUE
C
      DO 80 IA=1,L
      TSFA(IA)=1.0
C
      IF(IA .GT. 2)THEN
       DO 20 II=3,IA
       TSFA(IA-II+2)= TSFA(IA-II+2)*B + TSFA(IA-II+1)
       TSFAB(IA-II+2) = TSFA(IA-II+2)
   20  CONTINUE
      ENDIF
C
      IF(IA .GE. 2)THEN
       TSFA(1) = TSFA(1)*B
       TSFAB(1) = TSFA(1)
      ENDIF
C
      DO 70 IB=1,M
      IAB=IA+IB-1
      TSFAB(IAB)=1.
C
C     Calculate TSF.
C
      IF(IAB .GT. 1)THEN
       DO 50 II=2,IAB
       TSF(IADR+IAB-II+1) = TSFAB(IAB-II+2)*FAC(IAB-II+2)*C(IAB)
   50  CONTINUE
      ENDIF
      TSF(IADR) = TSFAB(1)*C(IAB)
C
C     Update TSFAB. This must be done after calculation of TSF.
C
      IF(IAB .GT. 1)THEN
       DO 60 II=2,IAB
       TSFAB(IAB-II+2)    = TSFAB(IAB-II+2)*A + TSFAB(IAB-II+1)
   60  CONTINUE
      ENDIF
      TSFAB(1)  = TSFAB(1)*A
C
      IADR=IADR+IAB
   70 CONTINUE
   80 CONTINUE
C
   90 CONTINUE
C
  100 CONTINUE
      RETURN
      END
