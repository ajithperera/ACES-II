      SUBROUTINE RCL_F2L2AB(T2,Z2,FMIAA,FMIBB,POP1,POP2,VRT1,VRT2,
     &                  DISSYT,DISSYZ,NUMSYT,NUMSYZ,
     &                  NFSIZA,NFSIZB,LISTT,LISTZ,
     &                  IRREP,IUHF,TMP,IMULT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DISSYT,DISSYZ,DIRPRD,POP1,POP2,VRT1,VRT2
      DIMENSION T2(DISSYT,NUMSYT),Z2(DISSYZ,NUMSYZ),FMIAA(NFSIZA),
     &          FMIBB(NFSIZB)
      DIMENSION TMP(1),POP1(8),POP2(8),VRT1(8),VRT2(8)
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
C
      DATA AZERO,ONE,ONEM,TWO /0.0D0,1.0D0,-1.0D0,2.0D0/
C
C   GET L2 AMPLITUDES
C
      CALL GETLST(T2,1,NUMSYT,2,IRREP,LISTT)

#ifdef _DEBUG_LVLM
      call checksum("L2ab :",t2,dissyt*numsyt)
CSSS      Write(6,*) " L2(Ab,Ij) increment at entering F1AB"
CSSS      call output(t2,1,dissyt,1,numsyt,dissyt,numsyt,1)

      call checksum("fmiaa:",fmiaa,NFSIZA)
      call checksum("fmibb:",fmibb,NFSIZA)
CSSS      call output(fmiaa,1,NFSIZA,1,1,NFSIZA,1,1)
#endif
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
       IRREPI=DIRPRD(IRREPJ,IRREP)
C
       NOCCI=POP1(IRREPI)
C
       IF(MIN(NOCCI,NOCCJ).NE.0) THEN
C
       CALL XGEMM('N','T',DISSYT*NOCCI,NOCCJ,NOCCJ,ONEM,T2(1,JOFF), 
     &            NOCCI*DISSYT,FMIBB(IOFF),NOCCJ,
     &            AZERO,Z2(1,JOFF),NOCCI*DISSYZ)
       ENDIF
C
       JOFF=JOFF+NOCCJ*NOCCI
       IOFF=IOFF+NOCCJ*NOCCJ
90    CONTINUE
C
      IF(IUHF.EQ.0) THEN

C IN RHF THIS IS SIMPLY A TRANSPOSITION

       CALL GETLST(T2,1,NUMSYZ,1,IRREP,LISTZ)

#ifdef _DEBUG_LVL0
      write(6,*) "@-rcc_fijab"
      call checksum("z2old:",t2,NUMSYZ*DISSYZ)
      call checksum("z2add:",z2,NUMSYZ*DISSYZ)
#endif
       CALL VADD(T2,T2,Z2,NUMSYZ*DISSYZ,ONE)
       CALL SYMTR1(IRREP,POP1,POP2,DISSYZ,Z2,TMP,TMP(1+DISSYZ),
     &             TMP(1+2*DISSYZ))
       CALL SYMTR3(IRREP,VRT1,VRT2,DISSYZ,NUMSYZ,Z2,TMP,
     &             TMP(1+NUMSYZ),TMP(1+2*NUMSYZ))

       IF (IMULT .EQ. 2) THEN
       CALL VADD(Z2,Z2,T2,NUMSYZ*DISSYZ,ONE)
       ELSE
       CALL VADD(Z2,Z2,T2,NUMSYZ*DISSYZ,ONE)
       ENDIF 

       CALL PUTLST(Z2,1,NUMSYZ,1,IRREP,LISTZ)
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
        IRREPJ=DIRPRD(IRREPI,IRREP)
C
        NOCCJ=POP2(IRREPJ)
C
        IF(MIN(NOCCI,NOCCJ).NE.0) THEN
C
        CALL XGEMM('N','T',DISSYT*NOCCJ,NOCCI,NOCCI,ONEM,T2(1,JOFF), 
     &             NOCCJ*DISSYT,FMIAA(IOFF),NOCCI,
     &             ONE,Z2(1,JOFF),NOCCJ*DISSYZ)
        ENDIF
C
        JOFF=JOFF+NOCCJ*NOCCI
        IOFF=IOFF+NOCCI*NOCCI
190    CONTINUE

       CALL SYMTR1(IRREP,POP2,POP1,DISSYZ,Z2,TMP,TMP(1+DISSYZ),
     &             TMP(1+2*DISSYZ))

       CALL GETLST(T2,1,NUMSYZ,1,IRREP,LISTZ)

#ifdef _DEBUG_LVL0
      write(6,*) "@-rcc_fij2ab"
      call checksum("z2old:",t2,NUMSYZ*DISSYZ)
      call checksum("z2add:",z2,NUMSYZ*DISSYZ)
#endif
       CALL VADD(Z2,Z2,T2,NUMSYZ*DISSYZ,ONE)
       CALL PUTLST(Z2,1,NUMSYZ,1,IRREP,LISTZ)

#ifdef _DEBUG_LVLM
      Write(6,*) " L2(Ab,Ij) increment at exit F1AB"
      call output(z2,1,dissyz,1,numsyz,dissyz,numsyz,1)
#endif

       ENDIF
       RETURN
       END
