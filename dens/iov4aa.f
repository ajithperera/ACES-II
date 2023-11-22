      SUBROUTINE IOV4AA(G,MAXSIZE,W,AIOV,FACT,ISPIN,POP,VRT,
     &                  DISSYT,NUMSYT,DISSYW,NUMSYW,LISTG,
     &                  LISTW,IRREP,TMP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYT,DISSYW,DISMAX,DISLEFT,DISREAD,DIRPRD,POP,VRT
      DIMENSION G(DISSYT,1),W(DISSYW,1),AIOV(1),
     &          POP(8),VRT(8),TMP(1),IPT(8)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      common /dropgeo/ ndrgeo
C
      DATA ONE /1.D0/
C
C  PICK UP THE INTEGRALS, WE NEED THEM ANYWAY
C
      CALL GETLST(W,1,NUMSYW,2,IRREP,LISTW)   
C
C EXPAND THE OCC. OCC. BLOCK SO THAT WE HAVE E<F,M,I
C
      CALL SYMEXP(IRREP,POP,DISSYW,W)
C
C  DECIDE ABOUT THE ALGORITHM
C
      IF(MAXSIZE.GE.DISSYT*NUMSYT) THEN
C
C  FULL IN CORE ALGORITHM
C
C      PICK UP THE G AMPLITUDES
C
       if (ndrgeo.eq.0) then
         CALL GETLST(G,1,NUMSYT,1,IRREP,LISTG)
       else  
         CALL GETGO4U(G,TMP,1,NUMSYT,1,IRREP,LISTG,ispin,listg,dissyt)
       endif 
C
C  TRANSPOSE THE LAST TWO INDICES (  E<F,A,M --> E<F,M,A)
C
       CALL SYMTR1(IRREP,VRT,POP,DISSYT,G,TMP,TMP(1+DISSYT),
     &             TMP(1+2*DISSYT))
C
C  PERFORM MULTIPLICATION 
C
C  JOFFG AND JOFFW OFFSET IN THE RIGHTMOST BLOCKS OF G AND W
C  IOFF OFFSET IN IOV
C
       IOFF=1
       JOFFG=1
       JOFFW=1
       DO 90 IRREPJ=1,NIRREP
C          
C GET NUMBER OF  OCCUPIED AND VIRTUAL ORBITALS FOR IRREPJ     
C
        NOCCJ=POP(IRREPJ)
        NVRTJ=VRT(IRREPJ)
C
C DETERMINE IRREPI WHOSE DIRECT PRODUCT WITH IRREPJ GIVES IRREP
C
        IRREPI=DIRPRD(IRREP,IRREPJ)
C
C GET NUMBER OF OCCUPIED ORBITALS FOR IRREPI
C
        NOCCI=POP(IRREPI)
C
C IF NOCCI, NOCCJ, OR NVRTJ EQUAL ZERO, NOTHING TO COMPUTE
C
        IF(MIN(NOCCI,NVRTJ,NOCCJ).NE.0) THEN
C
          CALL XGEMM('T','N',NVRTJ,NOCCJ,NOCCI*DISSYT,FACT,
     &              G(1,JOFFG),NOCCI*DISSYT,W(1,JOFFW),
     &              NOCCI*DISSYT,ONE,AIOV(IOFF),NVRTJ)
        ENDIF
C
C  UPDATE OFFSETS
C
        JOFFG=JOFFG+NOCCI*NVRTJ
        JOFFW=JOFFW+NOCCI*NOCCJ
        IOFF=IOFF+NVRTJ*NOCCJ
90     CONTINUE
C
      ELSE
C
C   WE HAVE TO IT OUT OF CORE
C
C  GET OFFSETS FOR AIOV
C
       IPT(1)=1
       DO 150 IRREPJ=1,NIRREP-1
        IPT(IRREPJ+1)=IPT(IRREPJ)+POP(IRREPJ)*VRT(IRREPJ)
150    CONTINUE
C   
C  DETERMINE MAXIMUM NUMBER OF DISTRIBUTIONS WHICH FIT INTO CORE
C
       MAXDIS=MAXSIZE/DISSYT
C
C SET OFFSETS FOR GAMMA LIST AND INTEGRALS
C
       IOFFSET=1
       JOFFW=1
       IF (NDRGEO.NE.0) THEN
        CALL GETINO4U(TMP,NUMSYT,DISSYT,LISTG,IRREP,ISPIN)
        ITIN = 1 + NUMSYT+1 
        NIKI = 1
        IODRSET = 1
       ENDIF 
C
C LOOP OVER THE IRREPS OF THE LAST INDEX
C
       DO 200 IRREPJ=1,NIRREP
C
        NOCCJ=POP(IRREPJ)
        IRREPI=DIRPRD(IRREP,IRREPJ)
        NOCCI=POP(IRREPI)
        NVRTI=VRT(IRREPI)
        IF(MIN(NOCCI,NVRTI,NOCCJ).NE.0) THEN
C
C DETERMINE MAXIMUM NUMBER OF (E<F,A) BLOCKS WHICH CAN BE HELD IN CORE
C
         DISMAX=MAXDIS/NVRTI 
         IF(DISMAX.LE.0) STOP 'IOV4AA'
C
C GET NUMBER OF (E<F,A) BLOCKS WHICH HAVE TO BE READ 
C
         DISLEFT=NOCCJ
C
10       CONTINUE
C
C DETERMINE NUMBER OF (E<F,A) BLOCKS WHICH ARE READ DURING THIS PASS
C
          DISREAD=MIN(DISLEFT,DISMAX)
          DISLEFT=DISLEFT-DISREAD
C
C GET THE DISTRIBUTIONS FROM DISK
C
          if (ndrgeo.eq.0) then 
           CALL GETLST(G,IOFFSET,DISREAD*NVRTI,1,IRREP,LISTG)
          else  
           CALL GETGO4UO(G,TMP(ITIN),IODRSET,DISREAD*NVRTI,1,IRREP,
     x              LISTG,ispin,listg,dissyt,TMP,NIKI)
          endif 
C
C UPDATE OFFSET
C
          IOFFSET=IOFFSET+DISREAD*NVRTI 
C
C  LOOP OVER ALL DISTRIBUTIONS AND PERFORM MULTIPLICATION
C
          IOFF=IPT(IRREPI)
          JOFFG=1
C
          DO 250 NUM=1,DISREAD
C
           CALL XGEMM('T','N',NVRTI,NOCCI,DISSYT,-FACT,
     &                G(1,JOFFG),DISSYW,W(1,JOFFW),DISSYT,
     &                ONE,AIOV(IOFF),NVRTI)
C
C  UPDATE OFFSETS
C
           JOFFW=JOFFW+NOCCI
           JOFFG=JOFFG+NVRTI
C 
250       CONTINUE
C
C  IF NOT ALL (E<F,A) BLOCKS HAVE BEEN PROCESSED, GO BACK TO 10
C
         IF(DISLEFT.NE.0) GO TO 10
C
        ELSE
C
C UPDATE OFFSETS IN THE CASE NOTHING HAS BEEN DONE
C
         JOFFW=JOFFW+NOCCI*NOCCJ
         IOFFSET=IOFFSET+NOCCJ*NVRTI
        ENDIF
C  
200    CONTINUE
C
      ENDIF
C
      RETURN
      END