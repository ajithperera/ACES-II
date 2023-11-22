      SUBROUTINE FMIAB1(T2,Z2,FMIAA,FMIBB,POP1,POP2,VRT1,VRT2,
     &                  DISSYT,DISSYZ,NUMSYT,NUMSYZ,
     &                  NFSIZA,NFSIZB,LISTT,LISTZ,
     &                  IRREP,IUHF,TMP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DISSYT,DISSYZ,DIRPRD,POP1,POP2,VRT1,VRT2
      DIMENSION T2(DISSYT,NUMSYT),Z2(DISSYZ,NUMSYZ),FMIAA(NFSIZA),
     &          FMIBB(NFSIZB)
      DIMENSION TMP(1),POP1(8),POP2(8),VRT1(8),VRT2(8)
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
C
      DATA AZERO,ONE,ONEM /0.0D0,1.0D0,-1.0D0/
C
C   GET T2 AND Z2 AMPLITUDES
C
      CALL GETLST(T2,1,NUMSYT,2,IRREP,LISTT)
      CALL ZERO(Z2,NUMSYZ*DISSYZ)
C
C   PERFORM MULTIPLICATION
C
      JOFF=1     
      IOFF=1
      DO 90 IRREPJ=1,NIRREP
C
       NOCCJ=POP2(IRREPJ)
C
       IF(NOCCJ.EQ.0) GO TO 90
C
       IRREPI=DIRPRD(IRREPJ,IRREP)
C
       NOCCI=POP1(IRREPI)
C
       IF(NOCCI.EQ.0) GO TO 80 
C
       CALL XGEMM('N','N',DISSYT*NOCCI,NOCCJ,NOCCJ,ONEM,T2(1,JOFF), 
     &            NOCCI*DISSYT,FMIBB(IOFF),NOCCJ,
     &            AZERO,Z2(1,JOFF),NOCCI*DISSYZ)
C
       JOFF=JOFF+NOCCJ*NOCCI
80     CONTINUE
       IOFF=IOFF+NOCCJ*NOCCJ
90    CONTINUE
C
      IF(IUHF.EQ.0) THEN
C
C   IN RHF THIS IS SIMPLY A TRANSPOSITION
C
       CALL GETLST(T2,1,NUMSYZ,1,IRREP,LISTZ)
       CALL SYMRHF(IRREP,VRT1,POP1,DISSYZ,Z2,TMP,TMP(1+DISSYZ),
     &             TMP(1+2*DISSYZ))
       CALL VADD(T2,T2,Z2,NUMSYZ*DISSYZ,ONE)
       CALL PUTLST(T2,1,NUMSYZ,1,IRREP,LISTZ)
CSSS       call checksum("fmi-ab", t2, NUMSYZ*DISSYZ)
    
      ELSE
       CALL SYMTR1(IRREP,POP1,POP2,DISSYT,T2,TMP,TMP(1+DISSYT),
     &             TMP(1+2*DISSYT))
       CALL SYMTR1(IRREP,POP1,POP2,DISSYZ,Z2,TMP,TMP(1+DISSYZ),
     &             TMP(1+2*DISSYZ))
       JOFF=1     
       IOFF=1
       DO 190 IRREPI=1,NIRREP
C
        NOCCI=POP1(IRREPI)
C
        IF(NOCCI.EQ.0) GO TO 190
C
        IRREPJ=DIRPRD(IRREPI,IRREP)
C
        NOCCJ=POP2(IRREPJ)
C
        IF(NOCCJ.EQ.0) GO TO 180 
C
        CALL XGEMM('N','N',DISSYT*NOCCJ,NOCCI,NOCCI,ONEM,T2(1,JOFF), 
     &             NOCCJ*DISSYT,FMIAA(IOFF),NOCCI,
     &             ONE,Z2(1,JOFF),NOCCJ*DISSYZ)
C
        JOFF=JOFF+NOCCJ*NOCCI
180    CONTINUE
       IOFF=IOFF+NOCCI*NOCCI
190    CONTINUE
       CALL SYMTR1(IRREP,POP2,POP1,DISSYZ,Z2,TMP,TMP(1+DISSYZ),
     &             TMP(1+2*DISSYZ))
       CALL GETLST(T2,1,NUMSYZ,1,IRREP,LISTZ)
       CALL VADD(Z2,Z2,T2,NUMSYZ*DISSYZ,ONE)
       CALL PUTLST(Z2,1,NUMSYZ,1,IRREP,LISTZ)
       ENDIF
       RETURN
       END
