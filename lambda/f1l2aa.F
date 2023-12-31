      SUBROUTINE F1L2AA(T2,Z2,T,Z,FEA,VRT,NVRTSQ,
     &                  DISSYT,DISSYZ,NUMSYT,NUMSYZ,
     &                  NFSIZ,LISTT,LISTZ,IRREP,TMP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DISSYT,DISSYZ,DIRPRD,VRT
      DIMENSION T2(DISSYT,NUMSYT),Z2(DISSYZ,NUMSYZ),FEA(NFSIZ),
     &          T(NUMSYT,NVRTSQ),Z(NUMSYZ,NVRTSQ)
      DIMENSION TMP(1),VRT(8)
C
      COMMON/SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &               DIRPRD(8,8)

      DATA AZERO,ONE,ONEM /0.0D0,1.0D0,-1.0D0/
C
C   GET L2 AMPLITUDES
C
      CALL GETLST(T2,1,NUMSYT,2,IRREP,LISTT)
      CALL TRANSP(T2,T,NUMSYT,DISSYT)
      CALL SYMEXP(IRREP,VRT,NUMSYT,T)
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
C
       JOFF=JOFF+NVRTJ*NVRTI
       IOFF=IOFF+NVRTJ*NVRTJ
90    CONTINUE
C
      CALL ASSYM(IRREP,VRT,NUMSYZ,NUMSYZ,T,Z)
      CALL TRANSP(T,T2,DISSYZ,NUMSYZ)
C
      CALL GETLST(Z2,1,NUMSYZ,1,IRREP,LISTZ)
      CALL VADD(Z2,Z2,T2,NUMSYZ*DISSYZ,ONE)
      CALL PUTLST(Z2,1,NUMSYZ,1,IRREP,LISTZ)
      RETURN
      END
