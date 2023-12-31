      SUBROUTINE IOV3BA(G,W,AIOV,FACT,ISPIN,POP1,POP2,VRT1,VRT2,
     &                  DISSYT,NUMSYT,DISSYW,NUMSYW,LISTG,LISTW,
     &                  IRREP,TMP,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYT,DISSYW,DIRPRD,POP1,POP2,VRT1,VRT2
      DIMENSION G(DISSYT,1),W(NUMSYW,1),AIOV(1),POP1(8),POP2(8),
     &          VRT1(8),VRT2(8),TMP(1)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      common /dropgeo/ ndrgeo
C
      DATA ONE /1.D0/
C
C      PICK UP THE INTEGRALS REQUIRED
C
      CALL GETLST(G,1,NUMSYW,2,IRREP,LISTW)
C
C  ISPIN=1   THE INTEGRALS ARE <Im//Ne> ORDERING (I,m,N,e)
C            THE REQUIRED ORDERING IS e,N,m,I
C            THEREFORE TRANSPOSE e and N, the IM AND Ne and
C            FINALLY I AND m
C  ISPIN=2   THE INTEGRALS ARE <Mi//En> ORDERING (M,i,E,n)
C            THE REQUIRED ORDERING IS E,n,M,i
C            TRANSPOSE ONLY En AND Mi
C
      IF(ISPIN.EQ.1) THEN
      CALL SYMTR1(IRREP,POP1,VRT2,DISSYW,G,TMP,TMP(1+DISSYW),
     &            TMP(1+2*DISSYW))
      ENDIF
C
      CALL TRANSP(G,W,NUMSYW,DISSYW)
C
      IF(ISPIN.EQ.1) THEN
      CALL SYMTR1(IRREP,POP1,POP2,NUMSYW,W,TMP,TMP(1+NUMSYW),
     &            TMP(1+2*NUMSYW))
      ENDIF
C
C  CHANGE ORDERING FROM eN,Am TO eN,mA OR En,aM TO En,Ma
C  
      if (ndrgeo.eq.0) CALL GETLST(G,1,NUMSYT,1,IRREP,LISTG)
      if (ndrgeo.eq.1) 
     x  CALL GETGO3Q(G,TMP,1,NUMSYT,1,IRREP,LISTG,ispin,listg,dissyt,1)
      CALL SYMTR1(IRREP,VRT1,POP2,DISSYT,G,TMP,TMP(1+DISSYT),
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
       NOCCJ=POP1(IRREPJ)
       NVRTJ=VRT1(IRREPJ)
C
C        DETERMINE IRREPI WHOSE DIRECT PRODUCT WITH IRREPJ GIVES IRREP
C
       IRREPI=DIRPRD(IRREP,IRREPJ)
C
C        GET NUMBER  OF OCCUPIED ORBITALS FOR IRREPI
C
       NOCCI=POP2(IRREPI)
C
C        IF NOCCI, NOCCJ, OR NVRTJ EQUAL ZERO, NOTHING TO COMPUTE
C
       IF(MIN(NOCCI,NVRTJ,NOCCJ).NE.0) THEN
C
          CALL XGEMM('T','N',NVRTJ,NOCCJ,NOCCI*DISSYT,FACT,
     &             G(1,JOFFG),NOCCI*DISSYT,W(1,JOFFW),
     &             NOCCI*DISSYT,ONE,AIOV(IOFF),NVRTJ)
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
