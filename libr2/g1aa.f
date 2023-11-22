      SUBROUTINE G1AA(L,T,T1,G,ISAME,ISPIN,POP,VRT,NVRTSQ,
     &                NVRT2SQ,DISSYW,DISSYT,NUMSYW,NUMSYT,
     &                LISTL,LISTT,IRREP,TMP,FACT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ISAME
      DOUBLE PRECISION L
      INTEGER DISSYT, DISSYW, DIRPRD,POP,VRT
      DIMENSION L(NUMSYW,NVRT2SQ),T(NUMSYT,NVRT2SQ)
      DIMENSION T1(DISSYT,NUMSYT),G(NVRTSQ)
      DIMENSION TMP(1)
      DIMENSION POP(8),VRT(8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
C
      DATA ONEM,ONE,HALF /-1.0D0,1.D0,0.5D0/
C
      FACTM=-FACT
C
C      PICK UP FIRST THE L AMPLITUDES
C      TRANSPOSE AND EXPAND THEM
C
      CALL GETLST(T1,1,NUMSYW,2,IRREP,LISTL)
      CALL TRANSP(T1,L,NUMSYW,DISSYW)
      CALL SYMEXP(IRREP,VRT,NUMSYW,L)
C 
C     FOR CC METHODS :
C     PICK UP THE T AMPLITUDES 
C     IF REQUIRED. TRANSPOSE AND EXPAND THE TAUS
C
      IF(.NOT.ISAME) THEN
       CALL GETLST(T1,1,NUMSYT,1,IRREP,LISTT)   
       CALL TRANSP(T1,T,NUMSYT,DISSYT)
       CALL SYMEXP(IRREP,VRT,NUMSYT,T)
      ELSE
C
C  FOR MBPT(4) SIMPLY COPY L TO T SINCE THEY ARE THE SAME  
C
c YAU : old
c      CALL ICOPY(NUMSYT*NVRT2SQ*IINTFP,L,1,T,1)
c YAU : new
       CALL DCOPY(NUMSYT*NVRT2SQ,L,1,T,1)
c YAU : end
      ENDIF
C 
      IOFF=1
      JOFF=1
      DO 90 IRREPJ=1,NIRREP
C          
C      GET NUMBER OF VIRTUAL ORBITALS FOR IRREPJ     
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
C        IF ZERO, NOTHING TO COMPUTE
C
       IF(MIN(NVRTJ,NVRTI).NE.0) THEN
C
        CALL XGEMM('T','N',NVRTJ,NVRTJ,NVRTI*NUMSYW,FACTM,
     &              T(1,JOFF),NVRTI*NUMSYW,L(1,JOFF),
     &              NVRTI*NUMSYW,ONE,G(IOFF),NVRTJ)
       ENDIF
C
       JOFF=JOFF+NVRTI*NVRTJ
       IOFF=IOFF+NVRTJ*NVRTJ
C
90    CONTINUE
C
      RETURN
      END