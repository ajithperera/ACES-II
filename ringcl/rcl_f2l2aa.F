      SUBROUTINE RCL_F2L2AA(T2,Z2,FMI,POP,VRT,NOCCSQ,DISSYT,DISSYZ,
     &                      NUMSYT,NUMSYZ,NFSIZ,LISTT,LISTZ,IRREP,
     &                      TMP,ISPIN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DISSYT,DISSYZ,DIRPRD,POP
      DIMENSION T2(DISSYT,NOCCSQ),Z2(DISSYZ,NOCCSQ),FMI(NFSIZ)
      DIMENSION TMP(1),POP(8),VRT(8)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
C
      DATA AZERO,ONE,ONEM,TWO/0.0D0,1.0D0,-1.0D0,2.0D0/
C
C   GET L2  AMPLITUDES
C
      CALL GETLST(T2,1,NUMSYT,2,IRREP,LISTT)

#ifdef _DEBUG_LVLM
      if (ispin .eq. 1) then 
CSSS      call checksum("L2 :",t2,dissyt*numsyt)
      Write(6,*) " L2(AB,IJ) increment at entering F1AA"
      call output(t2,1,dissyt,1,numsyt,dissyt,numsyt,1) 
CSSS      call checksum("fmi:",fmi,NOCCSQ)
CSSS      call output(fmi,1,NFSIZ,1,1,NFSIZ,1,1)
      endif 
#endif 
      IF (.NOT. (DRCCD .OR. RCCD)) THEN
         CALL SYMEXP(IRREP,POP,DISSYT,T2)
      ENDIF 
      CALL ZERO(Z2,DISSYZ*NOCCSQ)
C
C   PERFORM MULTIPLICATION
C
      JOFF=1     
      IOFF=1
      DO 90 IRREPJ=1,NIRREP
C
       NOCCJ=POP(IRREPJ)
C
       IRREPI=DIRPRD(IRREPJ,IRREP)
C
       NOCCI=POP(IRREPI)
C
       IF(MIN(NOCCI,NOCCJ).NE.0) THEN 
C
       CALL XGEMM('N','T',DISSYT*NOCCI,NOCCJ,NOCCJ,ONEM,T2(1,JOFF), 
     &            NOCCI*DISSYT,FMI(IOFF),NOCCJ,
     &            AZERO,Z2(1,JOFF),NOCCI*DISSYZ)
       ENDIF
C
       JOFF=JOFF+NOCCJ*NOCCI
       IOFF=IOFF+NOCCJ*NOCCJ
90    CONTINUE
C
      IF (.NOT. (DRCCD .OR. RCCD)) THEN
         CALL ASSYM(IRREP,POP,DISSYZ,DISSYZ,T2,Z2)
         CALL GETLST(Z2,1,NUMSYZ,1,IRREP,LISTZ)
         CALL VADD(Z2,Z2,T2,NUMSYZ*DISSYZ,ONE)
         CALL PUTLST(Z2,1,NUMSYZ,1,IRREP,LISTZ)
      ELSE
         CALL GETLST(T2,1,NUMSYZ,1,IRREP,LISTZ)
#ifdef _DEBUG_LVL0
      write(6,*) "@-rcc_fijaa"
      call checksum("z2old:",t2,NUMSYZ*DISSYZ)
      call checksum("z2add:",z2,NUMSYZ*DISSYZ)
#endif
         CALL VADD(T2,T2,Z2,NUMSYZ*DISSYZ,ONE)
         CALL SYMTR1(IRREP,POP,POP,DISSYZ,Z2,TMP,TMP(1+DISSYZ),
     &               TMP(1+2*DISSYZ))
         CALL SYMTR3(IRREP,VRT,VRT,DISSYZ,NUMSYZ,Z2,TMP,
     &               TMP(1+NUMSYZ),TMP(1+2*NUMSYZ))
         CALL VADD(Z2,Z2,T2,NUMSYZ*DISSYZ,ONE)
         CALL PUTLST(Z2,1,NUMSYZ,1,IRREP,LISTZ)

#ifdef _DEBUG_LVLM
      if (ispin .eq. 1) then
      Write(6,*) " L2(AB,IJ) increment at exit F1AA"
      call output(z2,1,dissyz,1,numsyz,dissyz,numsyz,1)
      endif
#endif
      ENDIF 

      RETURN
      END
