      SUBROUTINE GET_RNORM(IRREPX, NSIZEC, SR, TMP, MAXCOR, 
     &                     IUHF, RNORM)
C
C THIS ROUTINE PROPERLY NORMALIZES THE RIGHT-HAND EOM-CC EIGENVECTOR
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LAMTHERE, PRINT
      PARAMETER (TOL=1.D-5)
c
      DIMENSION SR(NSIZEC),TMP(MAXCOR)
      COMMON/LAMSTATE/LAMTHERE
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/FLAGS/IFLAGS(100)
C
      DATA ZILCH,ONE /0.0D0,1.0D0/
C
C THE FOLLOWING CODE ONLY RUNS FOR IRREPX=1
C
      IF (IRREPX .EQ. 1) THEN
C
C CALCULATE CONTRIBUTION FROM REFERENCE FUNCTION (S0)
C           -1 
C S0 = W   [ F(ai) C(ai) + W(ab,ij) C(ab,ij) ]
C
         CALL LOADVEC1(IRREPX, TMP, MAXCOR, IUHF, 93, 0, 13,
     &                 LENTOT, IUHF.EQ.0)
C
         ROOT = ONE
C
         RNORM0 = SDOT(LENTOT, TMP, 1, SR, 1)/ROOT
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6, '(1x,a,a,F10.5)') "The overlap with zeroth order ",
     &                        "operator:", RNORM0
#endif
C
C CALCULATE OVERLAP WITH LAMBDA VECTOR [SHOULD BE ZERO]
C
         CALL LOADVEC1(IRREPX, TMP, MAXCOR, IUHF, 190, 0, 143,
     &                 LENTOT, IUHF.EQ.0)
C
         RNORM = SDOT(LENTOT, TMP, 1, SR, 1)
C
      ENDIF
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6, '(1x,a,F10.5)') "The overlap with Lambda vector:", RNORM
      Write(6,*)
#endif
C
      RETURN
      END
