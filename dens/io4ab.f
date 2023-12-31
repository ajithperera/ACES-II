      SUBROUTINE IO4AB(G,W,AIOO,FACT,ISPIN,POP1,POP2,VRT1,VRT2,
     &                  DISSYT,NUMSYT,DISSYW,NUMSYW,LISTG,LISTW,IRREP,
     &                  TMP,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYT,DISSYW,DIRPRD,POP1,POP2,VRT1,VRT2
      DIMENSION G(DISSYT,1),W(DISSYW,1),AIOO(1000),TMP(1)
      DIMENSION POP1(8),POP2(8),VRT1(8),VRT2(8) 
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      common /dropgeo/ ndrgeo
C
      DATA ONE,ONEM,TWO /1.0D0,-1.D0,2.D0/
C
C      PICK UP FIRST THE G AMPLITUDES AND INTEGRALS
C
      if (ndrgeo.eq.0) CALL GETLST(G,1,NUMSYT,1,IRREP,LISTG)
      if (ndrgeo.eq.1) then 
       if (iuhf.eq.1 .and. ispin.eq.1) then
        CALL GETGO4R(G,W,1,NUMSYT,1,IRREP,LISTG,ISPIN,listw,dissyt)
       else 
        CALL GETGO4(G,W,1,NUMSYT,1,IRREP,LISTG,ISPIN,listw,dissyt)
       endif 
      endif 
C
      CALL GETLST(W,1,NUMSYW,2,IRREP,LISTW)
C
C  SPIN ADAPTED CODE FOR RHF
C
      IF(IUHF.EQ.0) THEN
       CALL SPINAD3(IRREP,VRT1,DISSYW,NUMSYW,W,TMP,TMP(1+NUMSYW))
      ENDIF
C
C  TRANSPOSE THE LAST TWO INDICES IN THE AA CASE (ONLY IF UHF)
C
      IF(ISPIN.EQ.1.AND.IUHF.EQ.1) THEN
        CALL SYMTR1(IRREP,POP1,VRT2,DISSYT,G,TMP,TMP(1+DISSYT),
     &              TMP(1+2*DISSYT))
        CALL SYMTR1(IRREP,POP1,VRT2,DISSYW,W,TMP,TMP(1+DISSYW),
     &              TMP(1+2*DISSYW))
      ENDIF
C
C  PERFORM MULTIPLICATION
C
C  JOFF OFFSET IN THE RIGHTMOST BLOCK OF G AND W
C  IOFF OFFSET IN AIVV
C
      IOFF=1
      JOFF=1
      DO 90 IRREPI=1,NIRREP
C          
C        GET NUMBER OF OCCUPIED ORBITALS FOR IRREPJ
C
       NOCCI=POP1(IRREPI)
C
C        DETERMINE IRREPI WHOSE DIRECT PRODUCT WITH IRREPJ GIVES IRREP
C
       IRREPJ=DIRPRD(IRREP,IRREPI)
C
C        GET NUMBER OF VIRTUAL ORBITALS FOR IRREPI
C
       NVRTJ=VRT2(IRREPJ)
C
C        IF NVRTJ OR NOCCI EQUAL ZERO, NOTHING TO COMPUTE
C
       IF(MIN(NOCCI,NVRTJ).NE.0) THEN
C
          CALL XGEMM('T','N',NOCCI,NOCCI,DISSYT*NVRTJ,FACT,
     &               W(1,JOFF),NVRTJ*DISSYT,G(1,JOFF),
     &               NVRTJ*DISSYT,ONE,AIOO(IOFF),NOCCI)
       ENDIF
C
C  UPDATE THE OFFSETS
C
       JOFF=JOFF+NOCCI*NVRTJ
       IOFF=IOFF+NOCCI*NOCCI
90    CONTINUE
C
      RETURN
      END
