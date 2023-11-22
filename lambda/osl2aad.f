      SUBROUTINE OSL2AAD(T,L,G,T1A,T1B,ISPIN,POPT1,POPT2,POPL1,POPL2,
     $   VRT1,VRT2,NOCCSQ,DISSYW,DISSYT,NUMSYW,NUMSYT,LISTL,IRREP,
     $   TMP,FACT,STRING)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION L
      CHARACTER*4 STRING
      INTEGER DISSYT, DISSYW, DIRPRD,POPT1,POPL1,POPL2,VRT1,VRT2,POPT2
      DIMENSION L(DISSYW,NUMSYW),T(DISSYT,NUMSYT),G(NOCCSQ)
      DIMENSION TMP(*)
      DIMENSION POPT1(8),POPL1(8),POPL2(8),VRT1(8),VRT2(8),POPT2(8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
C
      DATA ONE,ONEM,TWO /1.0D0,-1.D0,2.D0/,HALF /0.5D0/
C
      CALL FSGET(L,1,NUMSYW,1,IRREP,LISTL,STRING)
C
      CALL ZERO(T,NUMSYT*DISSYT)
      IF(ISPIN.EQ.1) THEN
         CALL FTAU(T,T1A,T1B,DISSYT,NUMSYT,POPT1,POPT2,VRT1,VRT2,IRREP,
     $        3,ONE)
      ELSE
         CALL FTAU(T,T1B,T1A,DISSYT,NUMSYT,POPT2,POPT1,VRT2,VRT1,IRREP,
     $        3,ONE)
      ENDIF
C
      IF(ISPIN.EQ.2) THEN 
         CALL SYMTR1(IRREP,POPL2,POPL1,DISSYW,L,TMP,TMP(1+DISSYW),
     &      TMP(1+2*DISSYW))
      ENDIF
C
      JOFFL=1
      JOFFT=1
      IOFF=1
      DO 90 IRREPI=1,NIRREP
C
         NOCCL=POPL2(IRREPI)
         NOCCT=POPT2(IRREPI)
C
         IRREPJ=DIRPRD(IRREP,IRREPI)
         NOCCJ=POPL1(IRREPJ)
C
         IF(MIN(NOCCL,NOCCT,NOCCJ).NE.0) THEN
C
            CALL XGEMM('T','N',NOCCT,NOCCL,DISSYW*NOCCJ,FACT*ONEM,
     &         T(1,JOFFT),NOCCJ*DISSYW,L(1,JOFFL),      
     &         NOCCJ*DISSYW,ONE,G(IOFF),NOCCT)
         ENDIF
C
         JOFFL=JOFFL+NOCCJ*NOCCL
         JOFFT=JOFFT+NOCCJ*NOCCT
         IOFF=IOFF+NOCCT*NOCCL
 90   CONTINUE
C
      RETURN
      END
