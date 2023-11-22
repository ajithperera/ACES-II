      SUBROUTINE G2AA(T,L,G,ISAME,ISPIN,POP,VRT,NOCCSQ,
     &                NOCC2SQ,DISSYW,DISSYT,NUMSYW,NUMSYT,LISTL,
     &                LISTT,IRREP,TMP,FACT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ISAME 
      DOUBLE PRECISION L 
      INTEGER DISSYT,DISSYW,DIRPRD,POP,VRT
      DIMENSION L(DISSYW,NOCC2SQ),T(DISSYT,NOCC2SQ),G(NOCCSQ)
      DIMENSION TMP(1),POP(8),VRT(8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
C
      DATA ONEM,ONE,HALF /-1.0D0,1.D0,0.5D0/
C
C      PICK UP FIRST THE L AMPLITUDES AND EXPAND THEM
C
      CALL GETLST(L,1,NUMSYW,2,IRREP,LISTL)
      CALL SYMEXP(IRREP,POP,DISSYW,L)
C
C  FOR MBPT4 COPY W TO T
C
      IF(ISAME) THEN
c YAU : old
c      CALL ICOPY(IINTFP*NOCC2SQ*DISSYT,L,1,T,1)
c YAU : new
       CALL DCOPY(NOCC2SQ*DISSYT,L,1,T,1)
c YAU : end
      ELSE
C
C  FOR CC GET T AMPLITUDES
C
      CALL GETLST(T,1,NUMSYT,1,IRREP,LISTT)
C
C  EXPAND T AMPLITUDES
C
      CALL SYMEXP(IRREP,POP,DISSYT,T)
      ENDIF
C
      JOFF=1
      IOFF=1
      DO 90 IRREPJ=1,NIRREP
C          
C        GET OCCUPATION NUMBER FOR JRREP     
C
       NOCCJ=POP(IRREPJ)
C
C        DETERMINE IRREPI WHOSE DIRECT PRODUCT WITH JRREP GIVES IRREP
C
       IRREPI=DIRPRD(IRREP,IRREPJ)
C
C        GET OCCUPATION NUMBER FOR IRREPI
C
       NOCCI=POP(IRREPI)
C
C        IF ZERO, NOTHING TO COMPUTE
C
       IF(MIN(NOCCJ,NOCCI).NE.0) THEN
C
        CALL XGEMM('T','N',NOCCJ,NOCCJ,DISSYW*NOCCI,FACT,
     &             T(1,JOFF),NOCCI*DISSYW,L(1,JOFF),      
     &             NOCCI*DISSYW,ONE,G(IOFF),NOCCJ)
       ENDIF
C
       JOFF=JOFF+NOCCJ*NOCCI
       IOFF=IOFF+NOCCJ*NOCCJ
90    CONTINUE
C
      RETURN
      END