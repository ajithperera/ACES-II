      SUBROUTINE DOOAA(T1,T2,DOO,FACT,ISPIN,POP,VRT,DISSYT,
     &                 NUMSYT,LISTT1,LISTT2,LISTT3,IRREP,IFLAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYT,DIRPRD,POP,VRT
C
      DIMENSION T1(DISSYT,NUMSYT),T2(DISSYT,NUMSYT),DOO(1),
     &          POP(8),VRT(8) 
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
C
      DATA ONEM,ONE,HALF,TWO /-1.0D0,1.D0,0.5D0,2.0D0/
C
C  PICK UP THE T1 AND T2 AMPLITUDES REQUIRED
C  CHECK IFLAG IN ORDER WHAT TO DO
C
      CALL GETLST(T1,1,NUMSYT,1,IRREP,LISTT1)
C
      IF(IFLAG.GT.1) THEN
C
C GENERAL CASE, GET T AND L AMPLITUDES (CC-METHODS)
C
       CALL GETLST(T2,1,NUMSYT,2,IRREP,LISTT2)

CSSS       Write(6,*) "The checksum of l2"
CSSS       call checksum("L2AAAA  ", T2, NUMSYT*DISSYT)
C
       IF(IFLAG.EQ.3) THEN
C
C FOR MBPT(3), FORM T[2] + 1/2 T[1], NOTE THAT THE VECTOR T2 CONTAINS
C T[1] + T[2], THEREFORE THE FACTOR IS ONEM
C
        CALL SSCAL(NUMSYT*DISSYT,TWO,T2,1)
        CALL SAXPY(NUMSYT*DISSYT,ONEM,T1,1,T2,1)
C
       ELSE IF(IFLAG.EQ.4) THEN
C
C FOR SECOND CALL MBPT(4), ...
C
        CALL SSCAL(NUMSYT*DISSYT,TWO,T1,1)
        CALL SAXPY(NUMSYT*DISSYT,ONE,T2,1,T1,1)
        CALL GETLST(T2,1,NUMSYT,2,IRREP,LISTT3)
C
       ELSE IF(IFLAG.EQ.5) THEN
C
C ROHF-MBPT(3), FORM  T[1] + T[2] + L[2]
C
        CALL SAXPY(NUMSYT*DISSYT,ONE,T1,1,T2,1)
        CALL GETLST(T1,1,NUMSYT,1,IRREP,LISTT3)
        CALL SAXPY(NUMSYT*DISSYT,ONE,T1,1,T2,1)
C
       ENDIF
      ELSE
C
C T1 AND T2 ARE IDENTICAL (MBPT(2) AND FIRST CALL MBPT(4)
C
c YAU : old
c      CALL ICOPY(NUMSYT*DISSYT*IINTFP,T1,1,T2,1)
c YAU : new
       CALL DCOPY(NUMSYT*DISSYT,T1,1,T2,1)
c YAU : end
C
      ENDIF
C
C EXPAND THE OCCUPIED-OCCUPIED BLOCK OF THE AMPLITUDES
C
      CALL SYMEXP(IRREP,POP,DISSYT,T1)
      CALL SYMEXP(IRREP,POP,DISSYT,T2)
C
C PERFORM MULTIPLICATION
C
C  JOFF OFFSET IN THE OCCUPIED-OCCUPIED BLOCK OF T1 AND T2
C  IOFF OFFSET IN DOO
C
      JOFF=1
      IOFF=1
      DO 90 IRREPJ=1,NIRREP
C          
C GET OCCUPATION NUMBER FOR JRREP     
C
       NOCCJ=POP(IRREPJ)
C
C DETERMINE IRREPI WHOSE DIRECT PRODUCT WITH JRREP GIVES IRREP
C
       IRREPI=DIRPRD(IRREP,IRREPJ)
C
C GET OCCUPATION NUMBER FOR IRREPI
C
       NOCCI=POP(IRREPI)
C
C IF NOCCI OR NOCCJ EQUAL ZERO, NOTHING TO COMPUTE
C
       IF(MIN(NOCCJ,NOCCI).NE.0) THEN
C
        CALL XGEMM('T','N',NOCCJ,NOCCJ,DISSYT*NOCCI,FACT,
     &              T1(1,JOFF),NOCCI*DISSYT,T2(1,JOFF),      
     &              NOCCI*DISSYT,ONE,DOO(IOFF),NOCCJ)
C
       ENDIF
C
C UPDATE THE OFFSETS
C
       JOFF=JOFF+NOCCJ*NOCCI
       IOFF=IOFF+NOCCJ*NOCCJ
90    CONTINUE
C
      RETURN
      END
