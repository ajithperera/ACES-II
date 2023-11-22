      SUBROUTINE DEWT1A(T1FIN,T1INI,DEN,EVAL,NT1,POP,VRT,NOCC,NVRT,
     1                  ISPIN,ITASK)
C
C     Subroutine for converting either a T1 like quantity to D1T1 or
C     D1T1 like quantity to T1. It must be called for each spin case.
C
C     ITASK = 0 => T1 --> D1T1. Otherwise, D1T1 --> T1.
C
C     This subroutine derives from a dependent of E4S.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
      DIMENSION T1FIN(NT1),EVAL(NOCC+NVRT),
     &          T1INI(NT1),DEN(NT1),POP(8),VRT(8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
C
C     GET ORBITAL ENERGIES FROM DISK
C
      IF(ISPIN.EQ.1)THEN
       CALL GETREC(20,'JOBARC','SCFEVALA',IINTFP*(NOCC + NVRT),
     &             EVAL)
      ELSE
       CALL GETREC(20,'JOBARC','SCFEVALB',IINTFP*(NOCC + NVRT),
     &             EVAL)
      ENDIF
C
C     LOOP OVER IRREPS 
C
      IND=0
      INDI=0
      INDA0=0
C
      DO  30 IRREP=1,NIRREP
C
      NOCCI=POP(IRREP)
      NVRTI=VRT(IRREP)
C
      DO  20 I=1,NOCCI
      INDI=INDI+1
      INDA=INDA0
      DO  10 IA=1,NVRTI
      INDA=INDA+1
      IND=IND+1
      DEN(IND) = EVAL(INDI) - EVAL(NOCC + INDA)
   10 CONTINUE
   20 CONTINUE
      INDA0=INDA0+NVRTI
   30 CONTINUE
C
C     --- T1 was passed, want D1T1 ---
C
      IF(ITASK.EQ.0)THEN
      DO 40 I=1,NT1
      T1FIN(I)=T1INI(I) * DEN(I)
   40 CONTINUE
      ENDIF
C
C     --- D1T1 was passed, want T1 ---
C
      IF(ITASK.NE.0)THEN
      DO 50 I=1,NT1
      T1FIN(I)=T1INI(I) / DEN(I)
   50 CONTINUE
      ENDIF
C
      RETURN
      END
