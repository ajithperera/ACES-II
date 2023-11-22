      SUBROUTINE IV4ABO(G,W,AIVV,FACT,ISPIN,POP1,POP2,VRT1,VRT2,
     &                  DISSYT,NUMSYT,DISSYW,NUMSYW,NDIS,LISTG,
     &                  LISTW,IRREP,TMP,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYT,DISSYW,DIRPRD,POP1,POP2,VRT1,VRT2
      DIMENSION G(DISSYT,NDIS),W(DISSYW,NUMSYW),AIVV(1000),TMP(1)
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
C PICK UP G-AMPLITUDES AS MEMORY PERMITS
C
      NLEFT=NUMSYW
      ISTART=1
      IF (NDRGEO.NE.0) THEN
        if (iuhf.eq.1 .and. ispin.eq.2) then 
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
       NREAD=MIN(NDIS,NLEFT)
       NLEFT=NLEFT-NREAD
       IEND=ISTART+NREAD-1
C
       if (ndrgeo.eq.0) then
         CALL GETLST(G,ISTART,NREAD,1,IRREP,LISTG)
       else 
         CALL GETGO4O(G,TMP(ITIN),ISTADR,NREAD,1,IRREP,LISTG,ispin,
     x                             listg,dissyt,TMP,NIKI) 
       endif 
C
C PERFORM MULTIPLICATION
C
C JOFF OFFSET IN THE RIGHTMOST BLOCK OF G AND W
C IOFF OFFSET IN AIVV
C
      IF(ISPIN.NE.1) THEN
C
C CODE FOR ISPIN EQUAL 2 (UHF ONLY)
C
      IOFF=1
      JOFF=1
C
      DO 190 IRREPI=1,NIRREP
C
C DETERMINE IRREPI WHOSE DIRECT PRODUCT WITH IRREPJ GIVES IRREP
C
       IRREPJ=DIRPRD(IRREP,IRREPI)
C
C GET NUMBER OF OCCUPIED ORBITALS FOR KRREP
C
       NOCCJ=POP2(IRREPJ)
C          
C GET NUMBER OF VIRTUAL ORBITALS FOR IRREPJ
C
       NVRTI=VRT1(IRREPI)
C
C IF NVRTI OR NOCCJ EQUAL ZERO, NOTHING TO COMPUTE
C
       IF(MIN(NOCCJ,NVRTI).NE.0) THEN
C
         DO 1200 IVRT1=1,NVRTI
          DO 1250 IVRT2=1,NVRTI
           DO 1300 IOCC=1,NOCCJ
            JOFF1=JOFF+IOCC+(IVRT1-1)*NOCCJ-1
            JOFF2=JOFF+IOCC+(IVRT2-1)*NOCCJ-1
            IOFF12=IOFF+IVRT1+(IVRT2-1)*NVRTI-1
            IF(JOFF2.GE.ISTART.AND.JOFF2.LE.IEND) THEN
             JOFF2=JOFF2-ISTART+1
              CALL XGEMM('T','N',1,1,DISSYT,FACT,
     &                  W(1,JOFF1),DISSYT,G(1,JOFF2),
     &                  DISSYT,ONE,AIVV(IOFF12),1)
            ENDIF
1300      CONTINUE
1250     CONTINUE
1200    CONTINUE
       ENDIF
C
C UPDATE OFFSETS
C
       IOFF=IOFF+NVRTI*NVRTI
       JOFF=JOFF+NVRTI*NOCCJ
C
190   CONTINUE
C
      ELSE
C
C CODE FOR ISPIN EQUAL 1 (RHF AND UHF)
C
      IOFF=1
      JOFF=1
C
      DO 90 IRREPJ=1,NIRREP
C
C DETERMINE IRREPI WHOSE DIRECT PRODUCT WITH IRREPJ GIVES IRREP
C
       IRREPI=DIRPRD(IRREP,IRREPJ)
C
       IOFF=1
       DO 91 IR=1,IRREPI-1
        IOFF=IOFF+VRT1(IR)*VRT1(IR)
91     CONTINUE
C          
C GET NUMBER OF VIRTUAL ORBITALS FOR IRREPJ
C
       NVRTI=VRT1(IRREPI)
C
C GET NUMBER OF OCCUPIED ORBITALS FOR KRREP
C
       NOCCJ=POP2(IRREPJ)
C
C IF NVRTI OR NOCCJ EQUAL ZERO, NOTHING TO COMPUTE
C
       IF(MIN(NOCCJ,NVRTI).NE.0) THEN
C
         DO 200 IVRT1=1,NVRTI
          DO 250 IVRT2=1,NVRTI
           DO 300 IOCC=1,NOCCJ
            JOFF1=JOFF+IVRT1+(IOCC-1)*NVRTI-1
            JOFF2=JOFF+IVRT2+(IOCC-1)*NVRTI-1
            IOFF12=IOFF+IVRT1+(IVRT2-1)*NVRTI-1
            IF(JOFF2.GE.ISTART.AND.JOFF2.LE.IEND) THEN
             JOFF2=JOFF2-ISTART+1
              CALL XGEMM('T','N',1,1,DISSYT,FACT,
     &                  W(1,JOFF1),DISSYT,G(1,JOFF2),
     &                  DISSYT,ONE,AIVV(IOFF12),NVRTI)
            ENDIF
300       CONTINUE
250      CONTINUE
200     CONTINUE
       ENDIF
C
C UPDATE OFFSETS
C
       JOFF=JOFF+NVRTI*NOCCJ
C
90    CONTINUE
C
      ENDIF
C
C UPDATE ISTART
C
      ISTART=ISTART+NREAD
C
      IF(NLEFT.NE.0) GO TO 1
C
      RETURN
      END
