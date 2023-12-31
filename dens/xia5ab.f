      SUBROUTINE XIA5AB(G,W,XIA,FACT,ISPIN,POP1,POP2,VRT1,VRT2,
     &                  DISSYT,NUMSYT,DISSYW,NUMSYW,LISTG,LISTW,
     &                  IRREP,TMP,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYT,DISSYW,DIRPRD,POP1,POP2,VRT1,VRT2
      DIMENSION G(DISSYT,NUMSYT),W(NUMSYW,DISSYW),XIA(1000),
     &          POP1(8),POP2(8),VRT1(8),VRT2(8),TMP(1)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      common /dropgeo/ ndrgeo
C
      DATA ONE /1.D0/
C
C ISPIN=1   THE INTEGRALS ARE <Ae//Nf> (ORDERING Ae,Nf)
C            CHANGE THE ORDERING TO  (fN,eA)
C            THIS INVOLVES A TRANSPOSITION OF N AND F
C            AND A TRANSPOSITION OF (Ae) AND (fN)
C            AND A TRANSPOSITION OF A AND e
C ISPIN=2   THE INTEGRALS ARE <Ea//Fn) (ORDERING Ea,Fn)
C            CHANGE THE ORDERING TO (Fn,Ea)
C            AND A TRANSPOSITION OF (Ea) and (Fn)
C
C RHF       IDENTICAL TO ISPIN=2
C
      if (ndrgeo.eq.0) then
        CALL GETLST(G,1,NUMSYT,1,IRREP,LISTG)
CSSS        print*,List_g
CSSS        call output(g,1,dissyt,1,numsyt,dissyt,numsyt,1)
      else 
        CALL GETGO3R(G,W,1,NUMSYT,1,IRREP,LISTG,ispin,listg,dissyt,1)
      endif 
C
C SPIN ADAPT FOR RHF
C
      IF(IUHF.EQ.0) THEN
       if (ndrgeo.eq.0) then
         CALL GETLST(W,1,NUMSYT,1,IRREP,123)
CSSS        print*,"G-123"
CSSS        call output(w,1,dissyt,1,numsyt,dissyt,numsyt,1)
       else  
         CALL GETGO3(W,TMP,1,NUMSYT,1,IRREP,123,ispin,123,dissyt)
       endif 
       CALL SAXPY(NUMSYT*DISSYT,ONE,W,1,G,1)
      ENDIF
C
C      PICK UP THE INTEGRALS REQUIRED
C
      CALL GETTRN(W,TMP,DISSYW,NUMSYW,1,IRREP,LISTW)
CSSS      print*,"w-30"
CSSS       call output(w,1,numsyw,1,dissyw,numsyw,dissyw,1)
C 
C REORDER INTEGRALS
C
      IF(ISPIN.EQ.1.AND.IUHF.EQ.1) THEN
       CALL SYMTR3(IRREP,POP1,VRT2,NUMSYW,DISSYW,W,TMP,
     &             TMP(1+DISSYW),TMP(1+2*DISSYW))
       CALL SYMTR1(IRREP,VRT1,VRT2,NUMSYW,W,TMP,TMP(1+NUMSYW),
     &             TMP(1+2*NUMSYW))
      ENDIF
C
C
C NOW ALL ARRAYS HAVE BEEN SET UP FOR THE MULTIPLICATION
C
C PERFORM MULTIPLICATION 
C
C JOFFG OFFSET IN THE OCCUPIED-VIRTUAL BLOCK OF G 
C JOFFW OFFSET IN THE VIRTUAL-VIRTUAL  BLOCK OF W
C IOFF OFFSET IN IOV
C
      IOFF=1
      JOFFG=1
      JOFFW=1
C
      DO 90 IRREPJ=1,NIRREP
C          
C GET NUMBER OF VIRTUAL ORBITALS FOR IRREPJ     
C
       NOCCJ=POP1(IRREPJ)
       NVRTJ=VRT1(IRREPJ)
C
C DETERMINE IRREPI WHOSE DIRECT PRODUCT WITH IRREPJ GIVES IRREP
C
       IRREPI=DIRPRD(IRREP,IRREPJ)
C
C GET NUMBER OF VIRTUAL ORBITALS FOR IRREPI
C
       NVRTI=VRT2(IRREPI)
C
C IF NVRTI, NOCCJ, OR NVRTJ EQUAL ZERO, NOTHING TO COMPUTE
C
       IF(MIN(NVRTI,NVRTJ,NOCCJ).NE.0) THEN
C
          CALL XGEMM('T','N',NVRTJ,NOCCJ,NVRTI*DISSYT,FACT,
     &             W(1,JOFFW),NVRTI*DISSYT,G(1,JOFFG),
     &             NVRTI*DISSYT,ONE,XIA(IOFF),NVRTJ)
       ENDIF
C
C UPDATE OFFSETS
C
       JOFFG=JOFFG+NVRTI*NOCCJ
       JOFFW=JOFFW+NVRTI*NVRTJ
       IOFF=IOFF+NVRTJ*NOCCJ
90    CONTINUE
C
      RETURN
      END
