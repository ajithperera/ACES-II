        SUBROUTINE DVVAA(T1,T2,T,DVV,FACT,ISPIN,POP,VRT,DISSYT,
     &                   NUMSYT,LISTT1,LISTT2,LISTT3,IRREP,IFLAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYT,DIRPRD,POP,VRT
C
      DIMENSION T(DISSYT,1)
      DIMENSION T1(NUMSYT,1),T2(NUMSYT,1),DVV(1),POP(8),VRT(8)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
C
      DATA ONE,ONEM,TWO /1.D0,-1.D0,2.D0/
C
C      PICK UP THE T1 AND T2 AMPLITUDES REQUIRED
C      CHECK IFLAG IN ORDER WHAT TO DO
C
      CALL GETLST(T,1,NUMSYT,1,IRREP,LISTT1)
      CALL TRANSP(T,T1,NUMSYT,DISSYT)
C
      IF(IFLAG.GT.1) THEN
C
C GENERAL CASE, GET T AND L AMPLITUDES (CC-METHODS)
C
       CALL GETLST(T,1,NUMSYT,2,IRREP,LISTT2)   
       CALL TRANSP(T,T2,NUMSYT,DISSYT)
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
        CALL GETLST(T,1,NUMSYT,2,IRREP,LISTT3)
        CALL TRANSP(T,T2,NUMSYT,DISSYT)
C
       ELSE IF(IFLAG.EQ.5) THEN
C
C ROHF-MBPT(3), FORM  T[1] + T[2] + L[2]
C
        CALL SAXPY(NUMSYT*DISSYT,ONE,T1,1,T2,1)
        CALL GETLST(T,1,NUMSYT,1,IRREP,LISTT3)
        CALL TRANSP(T,T1,NUMSYT,DISSYT)
        CALL SAXPY(NUMSYT*DISSYT,ONE,T1,1,T2,1)
C
       ELSE IF(IFLAG.EQ.6) THEN
C Regularized techniques, esp. for LinCC
C Scale T2 with reg/epsilon^2
        CALL GETLST(T,1,NUMSYT,2,IRREP,LISTT3)
        CALL TRANSP(T,T2,NUMSYT,DISSYT)
        call VECPRD(T1,T2,T1,NUMSYT*DISSYT)
        call getlst(T,1,NUMSYT,1,IRREP,LISTT2)
        CALL TRANSP(T,T2,NUMSYT,DISSYT)
       ENDIF
C
      ELSE
C
C T1 AND T2 ARE IDENTICAL (MBPT(2) AND FIRST CALL MBPT(4)
C
c YAU : old
c      CALL ICOPY(IINTFP*NUMSYT*DISSYT,T1,1,T2,1)
c YAU : new
       CALL DCOPY(NUMSYT*DISSYT,T1,1,T2,1)
c YAU : end
C
      ENDIF
C
      CALL SYMEXP(IRREP,VRT,NUMSYT,T1)
      CALL SYMEXP(IRREP,VRT,NUMSYT,T2)
C
C  PERFORM MULTIPLICATION 
C
C  JOFF OFFSET IN THE VIRTUAL-VIRTUAL BLOCK OF T1 AND T2
C  IOFF OFFSET IN DVV
C
      IOFF=1
      JOFF=1
      DO 90 IRREPJ=1,NIRREP
C          
C        GET NUMBER OF VIRTUAL ORBITALS FOR IRREPJ     
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
C        IF NVRTI OR NVRTJ EQUAL ZERO, NOTHING TO COMPUTE
C
       IF(MIN(NVRTI,NVRTJ).NE.0) THEN
C
        CALL XGEMM('T','N',NVRTJ,NVRTJ,NVRTI*NUMSYT,FACT,
     &             T1(1,JOFF),NVRTI*NUMSYT,T2(1,JOFF),
     &             NVRTI*NUMSYT,ONE,DVV(IOFF),NVRTJ)
       ENDIF
C
C  UPDATE OFFSETS
C
       JOFF=JOFF+NVRTI*NVRTJ
       IOFF=IOFF+NVRTJ*NVRTJ
90    CONTINUE
C
      RETURN
      END
