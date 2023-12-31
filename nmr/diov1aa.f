      SUBROUTINE DIOV1AA(G,W,T,DIOV,FACT,ISPIN,POP,VRT,
     &                   DISSYG,NUMSYG,DISSYW,NUMSYW,
     &                   LISTG,LISTW,IRREPGL,IRREPL,
     &                   IRREPR,TMP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYG,DISSYW,DIRPRD,POP,VRT
      DIMENSION T(1),G(NUMSYG,1),W(DISSYW,1),DIOV(1),
     &          POP(8),VRT(8),TMP(1)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
      DATA ONE /1.D0/
C
C PICK UP THE G AMPLITUDES AND INTEGRALS REQUIRED
C
      CALL GETLST(T,1,NUMSYG,1,IRREPGL,LISTG)
      CALL TRANSP(T,G,NUMSYG,DISSYG)
      CALL SYMEXP(IRREPL,VRT,NUMSYG,G)
      CALL GETLST(W,1,NUMSYW,2,IRREPR,LISTW)   
      CALL SYMTR1(IRREPR,POP,VRT,DISSYW,W,TMP,TMP(1+DISSYW),
     &            TMP(1+2*DISSYW))
C
C PERFORM MULTIPLICATION 
C
C JOFFG AND JOFFW OFFSET IN THE VIRTUAL-VIRTUAL BLOCK OF G AND W
C IOFF OFFSET IN IOV
C
      IOFF=1
      JOFFW=1
      DO 90 IRREPJR=1,NIRREP
C          
C GET NUMBER OF VIRTUAL ORBITALS FOR IRREPJ     
C
       NOCCJR=POP(IRREPJR)
C
C DETERMINE IRREPI WHOSE DIRECT PRODUCT WITH IRREPJ GIVES IRREP
C
       IRREPI=DIRPRD(IRREPR,IRREPJR)
C
C GET NUMBER  OF VIRTUAL ORBITALS FOR IRREPI
C
       NVRTI=VRT(IRREPI)
C
       IRREPJL=DIRPRD(IRREPI,IRREPL)
C
       NVRTJL=VRT(IRREPJL)
C
       JOFFG=1
       DO 89 IRREP=1,IRREPJL-1
        JOFFG=JOFFG+VRT(IRREP)*VRT(DIRPRD(IRREP,IRREPL))
89     CONTINUE
C
C IF NVRTI OR NOCCJR OR NVRTJL EQUAL ZERO, NOTHING TO COMPUTE
C
       IF(MIN(NVRTI,NVRTJL,NOCCJR).NE.0) THEN
C
        CALL XGEMM('T','N',NVRTJL,NOCCJR,NVRTI*NUMSYG,FACT,
     &             G(1,JOFFG),NVRTI*NUMSYG,W(1,JOFFW),
     &             NVRTI*NUMSYG,ONE,DIOV(IOFF),NVRTJL)
       ENDIF
C
C UPDATE OFFSETS
C
       JOFFW=JOFFW+NVRTI*NOCCJR
       IOFF=IOFF+NVRTJL*NOCCJR
C
90    CONTINUE
C
      RETURN
      END
