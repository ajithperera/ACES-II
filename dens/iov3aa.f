      SUBROUTINE IOV3AA(G,W,AIOV,FACT,ISPIN,POP,VRT,DISSYT,
     &                  NUMSYT,DISSYW,NUMSYW,LISTG,LISTW,IRREP,TMP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYT,DISSYW,DIRPRD,POP,VRT
      DIMENSION G(DISSYT,1),W(NUMSYW,1),AIOV(1),POP(8),VRT(8),TMP(1)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      common /dropgeo/ ndrgeo
C
      DATA ONE /1.D0/
C
C      PICK UP THE INTEGRALS REQUIRED
C
C NOTE THE INTEGRALS HAVE THE ORDERING : ( M < I ; N   E)
C THE REQUIRED ORDERING IS :             ( E   N ; M   I)
C THUS. THE FOLLWONG OPERATIONS ARE REQUIRED :
C
C  TRANSPOSE E AND N          ( ---> SYMTR1 CALL)
C  TRANSPOSE M,I WITH E,N     ( ---> TRANSP CALL)
C  EXPAND M < I TO FULL LIST  ( ---> SYMEXP CALL)
C
      CALL GETLST(G,1,NUMSYW,2,IRREP,LISTW)
      CALL SYMTR1(IRREP,POP,VRT,DISSYW,G,TMP,TMP(1+DISSYW),
     &            TMP(1+2*DISSYW))
      CALL TRANSP(G,W,NUMSYW,DISSYW)
      CALL SYMEXP(IRREP,POP,NUMSYW,W) 
C
C  PICK UP THE GAMMA AMPLITUDES
C
C  THEIR ORDERING IS AFTER THE RESORT E,N ; A,M
C  TRANSPOSE HERE THE LAST TWO INDICES ( ---> SYMTR1 CALL)
C
      if (ndrgeo.eq.0) CALL GETLST(G,1,NUMSYT,1,IRREP,LISTG)
      if (ndrgeo.eq.1)
     x     CALL GETGO3(G,TMP,1,NUMSYT,1,IRREP,LISTG,ispin,listg,dissyt)
      CALL SYMTR1(IRREP,VRT,POP,DISSYT,G,TMP,TMP(1+DISSYT),
     &            TMP(1+2*DISSYT))
C
C  NOW ALL ARRAYS HAVE BEEN SET UP FOR THE MULTIPLICATION
C
C  PERFORM MULTIPLICATION 
C
C  JOFFG OFFSET IN THE OCCUPIED-VIRTUAL BLOCK OF G 
C  JOFFW OFFSET IN THE OCCUPIED-OCCUPIED BLOCK OF W
C  IOFF OFFSET IN IOV
C
      IOFF=1
      JOFFG=1
      JOFFW=1
      DO 90 IRREPJ=1,NIRREP
C          
C        GET NUMBER OF VIRTUAL ORBITALS FOR IRREPJ     
C
       NOCCJ=POP(IRREPJ)
       NVRTJ=VRT(IRREPJ)
C
C        DETERMINE IRREPI WHOSE DIRECT PRODUCT WITH IRREPJ GIVES IRREP
C
       IRREPI=DIRPRD(IRREP,IRREPJ)
C
C        GET NUMBER  OF OCCUPIED ORBITALS FOR IRREPI
C
       NOCCI=POP(IRREPI)
C
C        IF NOCCI, NOCCJ OR NVRTJ EQUAL ZERO, NOTHING TO COMPUTE
C
       IF(MIN(NOCCI,NVRTJ,NOCCJ).NE.0) THEN
C
          CALL XGEMM('T','N',NVRTJ,NOCCJ,NOCCI*NUMSYT,FACT,
     &             G(1,JOFFG),NOCCI*NUMSYT,W(1,JOFFW),
     &             NOCCI*NUMSYT,ONE,AIOV(IOFF),NVRTJ)
       ENDIF
C
C  UPDATE OFFSETS
C
       JOFFG=JOFFG+NOCCI*NVRTJ
       JOFFW=JOFFW+NOCCI*NOCCJ
       IOFF=IOFF+NVRTJ*NOCCJ
90    CONTINUE
C
      RETURN
      END
