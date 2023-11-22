      SUBROUTINE IO4ABO(G,W,AIOO,FACT,ISPIN,POP1,POP2,VRT1,VRT2,
     &                  DISSYT,NUMSYT,DISSYW,NUMSYW,NDIS,
     &                  LISTG,LISTW,IRREP,TMP,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYT,DISSYW,DIRPRD,POP1,POP2,VRT1,VRT2
      DIMENSION G(DISSYT,NDIS),W(DISSYW,NUMSYW),AIOO(1000),TMP(1)
      DIMENSION POP1(8),POP2(8),VRT1(8),VRT2(8) 
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD  
      common /dropgeo/ ndrgeo
C
      DATA ONE,ONEM,TWO /1.0D0,-1.D0,2.D0/
C
C PICK UP INTEGRALS
C
      CALL GETLST(W,1,NUMSYW,2,IRREP,LISTW)
C
C  SPIN ADAPTED CODE FOR RHF
C
      IF(IUHF.EQ.0) THEN
C
       CALL SPINAD3(IRREP,VRT1,DISSYW,NUMSYW,W,TMP,TMP(1+NUMSYW))
C
      ENDIF
C
C PICK UP G-AMPLITUDES AS CORE MEMORY PERMITS 
C
      ISTART=1
      NLEFT=NUMSYW
      IF (NDRGEO.NE.0) THEN
        if (iuhf.eq.1 .and. ispin.eq.1) then 
         CALL GETINO4R(TMP,NUMSYT,DISSYT,LISTG,IRREP)
        else 
         CALL GETINO4 (TMP,NUMSYT,DISSYT,LISTG,IRREP)
        endif 
        ITIN = 1 + NUMSYT+1
        NIKI = 1
        ISTADR = 1
      ENDIF 
C      
1     CONTINUE 
C
C DETERMINE NUMBER OF DISTRIBUTIONS TO READ FROM DISK IN THIS PASS
C
       NREAD=MIN(NLEFT,NDIS)
       NLEFT=NLEFT-NREAD
       IEND=ISTART+NREAD-1 
C 
C READ DISTRIBUTIONS FROM DISK   
C 
      if (ndrgeo.eq.0) then
        CALL GETLST(G,ISTART,NREAD,1,IRREP,LISTG)
      else 
        CALL GETGO4O(G,TMP(ITIN),ISTADR,NREAD,1,IRREP,LISTG,ISPIN,
     x                listg,dissyt,TMP,NIKI)
      endif 
C
      IF(ISPIN.EQ.1.AND.IUHF.EQ.1) THEN
C
C  CODE FOR UHF IN CASE OF ISPIN EQUAL 1
C
C  PERFORM MULTIPLICATION
C
C  JOFF OFFSET IN THE RIGHTMOST BLOCK OF G AND W
C  IOFF OFFSET IN AIVV
C
      IOFF=1
      JOFF=1
C
      DO 190 IRREPJ=1,NIRREP
C
C GET NUMBER OF VIRTUAL ORBITALS FOR IRREPI
C
       NVRTJ=VRT2(IRREPJ)
C
C DETERMINE IRREPI WHOSE DIRECT PRODUCT WITH IRREPJ GIVES IRREP
C
       IRREPI=DIRPRD(IRREP,IRREPJ)
C
       IOFF=1
       DO 191 IR=1,IRREPI-1
        IOFF=IOFF+POP1(IR)*POP1(IR)
191    CONTINUE
C          
C GET NUMBER OF OCCUPIED ORBITALS FOR IRREPJ
C
       NOCCI=POP1(IRREPI)
C
C IF NVRTJ OR NOCCI EQUAL ZERO, NOTHING TO COMPUTE
C
       IF(MIN(NOCCI,NVRTJ).NE.0) THEN
C
        DO 1200 IOCC1=1,NOCCI
C
         DO 1250  IOCC2=1,NOCCI
          DO 1300 IVRT=1,NVRTJ
C
            JOFF1=JOFF+(IVRT-1)*NOCCI+IOCC1-1
            JOFF2=JOFF+(IVRT-1)*NOCCI+IOCC2-1
            IOFF12=IOFF+IOCC1+(IOCC2-1)*NOCCI-1
            IF(JOFF2.GE.ISTART.AND.JOFF2.LE.IEND) THEN
             JOFF2=JOFF2-ISTART+1
              CALL XGEMM('T','N',1,1,DISSYT,FACT,
     &                  W(1,JOFF1),DISSYT,G(1,JOFF2),
     &                  DISSYT,ONE,AIOO(IOFF12),1)
C
            ENDIF
1300      CONTINUE
1250     CONTINUE
1200    CONTINUE
       ENDIF
C
C UPDATE OFFSET
C
       JOFF=JOFF+NOCCI*NVRTJ
C
190   CONTINUE
C

      ELSE
C
C CODE FOR RHF AND UHF IN CASE OF ISPIN EQUAL 2
C
C PERFORM MULTIPLICATION
C
C  JOFF OFFSET IN THE RIGHTMOST BLOCK OF G AND W
C  IOFF OFFSET IN AIVV
C
      IOFF=1
      JOFF=1
C
      DO 90 IRREPI=1,NIRREP
C          
C GET NUMBER OF OCCUPIED ORBITALS FOR IRREPJ
C
       NOCCI=POP1(IRREPI)
C
C DETERMINE IRREPI WHOSE DIRECT PRODUCT WITH IRREPJ GIVES IRREP
C
       IRREPJ=DIRPRD(IRREP,IRREPI)
C
C GET NUMBER OF VIRTUAL ORBITALS FOR IRREPI
C
       NVRTJ=VRT2(IRREPJ)
C
C IF NVRTJ OR NOCCI EQUAL ZERO, NOTHING TO COMPUTE
C
       IF(MIN(NOCCI,NVRTJ).NE.0) THEN
C
        DO 200 IOCC1=1,NOCCI
C
         DO 250  IOCC2=1,NOCCI
          DO 300 IVRT=1,NVRTJ
C
            JOFF1=JOFF+(IOCC1-1)*NVRTJ+IVRT-1
            JOFF2=JOFF+(IOCC2-1)*NVRTJ+IVRT-1
            IOFF12=IOFF+IOCC1+(IOCC2-1)*NOCCI-1
            IF(JOFF2.GE.ISTART.AND.JOFF2.LE.IEND) THEN
             JOFF2=JOFF2-ISTART+1
             CALL XGEMM('T','N',1,1,DISSYT,FACT,
     &                  W(1,JOFF1),DISSYT,G(1,JOFF2),
     &                  DISSYT,ONE,AIOO(IOFF12),1)
            ENDIF
300       CONTINUE
250      CONTINUE
200     CONTINUE
       ENDIF
C
C  UPDATE THE OFFSETS
C
       JOFF=JOFF+NOCCI*NVRTJ
       IOFF=IOFF+NOCCI*NOCCI
C
90    CONTINUE
C
      ENDIF
C
C UPDATE ISTART
C
      ISTART=ISTART+NREAD
C
C IF THERE ARE DISTRIBUTIONS LEFT ON DISK, GO TO 1
C
      IF(NLEFT.NE.0) GO TO 1
C
      RETURN
      END
