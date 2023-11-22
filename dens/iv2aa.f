      SUBROUTINE IV2AA(G,W,T,MAXDIS,AIVV,FACT,ISPIN,POP,VRT,
     &                 DISSYT,NUMSYT,LISTG,LISTW,IRREP,TMP)
C
C     THIS VERSION HANDLES IN CORE AND OUT CORE LOGIC FOR
C     THE <AE//FG> G(BE,FG) CONTRACTION
C     
C  CODED AUGUST/90  JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYT,DIRPRD,POP,VRT
      DIMENSION T(1),G(1),W(1),AIVV(1),
     &          POP(8),VRT(8),TMP(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD 
      common /dropgeo/ ndrgeo
C
      DATA ONEM,ONE,HALF /-1.0D0,1.D0,0.5D0/
C
C      PICK UP THE G AMPLITUDES AND INTEGRALS REQUIRED
C
      IOFFSET=1
      NTDIS=NUMSYT
      IF (NDRGEO.NE.0) THEN
       CALL GETINV2U(TMP,NUMSYT,DISSYT,LISTG,IRREP,ISPIN)
       NIKI = 1
       IODRSET = 1
      ENDIF 
1     NRDIS=MIN(MAXDIS,NTDIS)
      NTDIS=NTDIS-NRDIS
C
      if (ndrgeo.eq.0) THEN
       CALL GETLST(T,IOFFSET,NRDIS,1,IRREP,LISTG)
      else  
       CALL GETGV2UO(T,G,IODRSET,NRDIS,1,IRREP,LISTG,ispin,listg,dissyt,
     x                TMP,NIKI)
      endif 
      CALL TRANSP(T,G,NRDIS,DISSYT)
      CALL SYMEXP(IRREP,VRT,NRDIS,G)
C
      CALL GETLST(T,IOFFSET,NRDIS,2,IRREP,LISTW)   
      CALL TRANSP(T,W,NRDIS,DISSYT)
      CALL SYMEXP(IRREP,VRT,NRDIS,W)
C
C UPDATE IOFFSET
C
      IOFFSET=IOFFSET+NRDIS
C
C  PERFORM MULTIPLICATION 
C
C  JOFF OFFSET IN THE VIRTUAL-VIRTUAL BLOCK OF G AND W
C  IOFF OFFSET IN DVV
C
      IOFF=1
      JOFF=1
      DO 90 IRREPJ=1,NIRREP
C          
C        GET NUMBER OF VIRTUAL ORBITALS FOR IRREPJ     
C
       NVRTJ=VRT(IRREPJ)
C
C        DETERMINE IRREPI WHOSE DIRECT PRODUCT WITH IRREPJ GIVES IRREP
C
       IRREPI=DIRPRD(IRREP,IRREPJ)
C
C        GET NUMBER  OF VIRTUAL ORBITALS FOR IRREPI
C
       NVRTI=VRT(IRREPI)
C
C        IF NVRTI OR NVRTJ EQUAL ZERO, NOTHING TO COMPUTE
C
       IF(MIN(NVRTI,NVRTJ).NE.0) THEN
C
        CALL XGEMM('T','N',NVRTJ,NVRTJ,NVRTI*NRDIS,FACT,
     &             W(JOFF),NVRTI*NRDIS,G(JOFF),
     &             NVRTI*NRDIS,ONE,AIVV(IOFF),NVRTJ)
       ENDIF
C
C  UPDATE OFFSETS
C
       JOFF=JOFF+NVRTI*NVRTJ*NRDIS
       IOFF=IOFF+NVRTJ*NVRTJ
90    CONTINUE
      IF(NTDIS.NE.0) GO TO 1
C
      RETURN
      END