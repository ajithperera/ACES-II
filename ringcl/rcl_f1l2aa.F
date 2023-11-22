      SUBROUTINE RCL_F1L2AA(T,T1,Z,FEA,VRT,POP,NVRTSQ,
     &                      DISSYT,DISSYZ,NUMSYT,NUMSYZ,
     &                      NFSIZ,LISTT,LISTZ,IRREP,TMP,ISPIN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DISSYT,DISSYZ,DIRPRD,VRT
      DIMENSION T1(DISSYT,NUMSYT),FEA(NFSIZ),
     &          T(NUMSYT,NVRTSQ),Z(NUMSYZ,NVRTSQ)
      DIMENSION TMP(1),VRT(8),POP(8)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
      COMMON/SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &               DIRPRD(8,8)
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
      DATA AZERO,ONE,ONEM,TWO/0.0D0,1.0D0,-1.0D0,2.0D0/
C
C   GET L2 AMPLITUDES
C
      CALL GETLST(T1,1,NUMSYT,2,IRREP,LISTT)
#ifdef _DEBUG_LVLM
      if (ispin .eq.1) then 
CSSS      call checksum("L2aa :",t1,dissyt*numsyt)
      Write(6,*) " L2(AB,IJ) increment at entering F1AA"
      call output(t1,1,dissyt,1,numsyt,dissyt,numsyt,1) 
CSSS      call checksum("fea:",fea,NVRTSQ)
      write(6,*) " Fea"
      call output(fea,1,NFSIZ,1,1,NFSIZ,1,1)
      endif 
#endif

      CALL TRANSP(T1,T,NUMSYT,DISSYT)
      IF (.NOT. (DRCCD .OR. RCCD)) THEN
         CALL SYMEXP(IRREP,VRT,NUMSYT,T)
      ENDIF 
      CALL ZERO(Z,NUMSYZ*NVRTSQ)
C
C      PERFORM MULTIPLICATION
C
      JOFF=1
      IOFF=1
      DO 90 IRREPJ=1,NIRREP
C
       NVRTJ=VRT(IRREPJ)
C
       IRREPI=DIRPRD(IRREPJ,IRREP)
C
       NVRTI=VRT(IRREPI)
C
       IF(MIN(NVRTI,NVRTJ).NE.0) THEN
C 
       CALL XGEMM('N','T',NUMSYT*NVRTI,NVRTJ,NVRTJ,ONE,T(1,JOFF),
     &            NVRTI*NUMSYT,FEA(IOFF),NVRTJ,
     &            AZERO,Z(1,JOFF),NVRTI*NUMSYZ)
       ENDIF
       JOFF=JOFF+NVRTJ*NVRTI
       IOFF=IOFF+NVRTJ*NVRTJ
90    CONTINUE
C
      IF (.NOT. (DRCCD .OR. RCCD)) THEN
         CALL ASSYM(IRREP,VRT,NUMSYZ,NUMSYZ,T1,Z)
         CALL TRANSP(T1,T,DISSYZ,NUMSYZ)
         CALL GETLST(Z,1,NUMSYZ,1,IRREP,LISTZ)
         CALL VADD(Z,Z,T,NUMSYZ*DISSYZ,ONE)
         CALL PUTLST(Z,1,NUMSYZ,1,IRREP,LISTZ)
      ELSE
        CALL TRANSP(Z,T1,DISSYZ,NUMSYZ)
        CALL GETLST(T,1,NUMSYZ,1,IRREP,LISTZ)
#ifdef _DEBUG_LVL0
        write(6,*) "@-rcl_fabaa"
        call checksum("z2old:",t,NUMSYZ*DISSYZ)
        call checksum("z2add:",t1,NUMSYZ*DISSYZ)
#endif
        CALL VADD(T,T,T1,NUMSYZ*DISSYZ,ONE)
        CALL SYMTR1(IRREP,POP,POP,DISSYZ,T1,TMP,TMP(1+DISSYZ),
     &             TMP(1+2*DISSYZ))
        CALL SYMTR3(IRREP,VRT,VRT,DISSYZ,NUMSYZ,T1,TMP,
     &             TMP(1+NUMSYZ),TMP(1+2*NUMSYZ))
        CALL VADD(T1,T1,T,NUMSYZ*DISSYZ,ONE)
        CALL PUTLST(T1,1,NUMSYZ,1,IRREP,LISTZ)

#ifdef _DEBUG_LVLM
      if (ispin .eq. 1) then
      Write(6,*) " L2(AB,IJ) increment at exit F1AA"
      call output(t1,1,dissyz,1,numsyz,dissyz,numsyz,1)
      endif
#endif
      ENDIF 

      RETURN
      END
