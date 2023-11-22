        SUBROUTINE DDVVAA(T1,T2,DDVV,FACT,ISPIN,POP,VRT,DISSYL,
     &                    NUMSYL,DISSYR,NUMSYR,LISTT1,LISTT2,
     &                    LISTT3,IRREPL,IRREPR,IFLAG,TMP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYL,DISSYR,DIRPRD,POP,VRT
      DIMENSION TMP(1)
      DIMENSION T1(NUMSYL,1),T2(NUMSYR,1),DDVV(1),POP(8),VRT(8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
C
      DATA ONE,ONEM,TWO /1.D0,-1.D0,2.D0/
C
C      PICK UP THE T1 AND T2 AMPLITUDES REQUIRED
C      CHECK ISAME IF THEY ARE THE SAME OR NOT
C
      CALL GETTRN(T1,TMP,DISSYL,NUMSYL,1,IRREPR,LISTT1)
      CALL GETTRN(T2,TMP,DISSYR,NUMSYR,1,IRREPR,LISTT2)
c       IF(IFLAG.EQ.3) THEN
c        CALL SSCAL(NUMSYT*DISSYT,TWO,T2,1)
c        CALL SAXPY(NUMSYT*DISSYT,ONEM,T1,1,T2,1)
c       ELSE IF(IFLAG.EQ.4) THEN
c        CALL SSCAL(NUMSYT*DISSYT,TWO,T1,1)
c        CALL SAXPY(NUMSYT*DISSYT,ONE,T2,1,T1,1)
c        CALL GETLST(T,1,NUMSYT,2,IRREP,LISTT3)
c        CALL TRANSP(T,T2,NUMSYT,DISSYT)
c       ENDIF
c      ELSE
c       CALL ICOPY(IINTFP*NUMSYT*DISSYT,T1,1,T2,1)
c      ENDIF
      CALL SYMEXP(IRREPL,VRT,NUMSYL,T1)
      CALL SYMEXP(IRREPR,VRT,NUMSYR,T2)
C
C  PERFORM MULTIPLICATION 
C
C  JOFF OFFSET IN THE VIRTUAL-VIRTUAL BLOCK OF T1 AND T2
C  IOFF OFFSET IN DVV
C
      IOFF=1
      JOFF2=1
C
      DO 90 IRREPJR=1,NIRREP
C          
C        GET NUMBER OF VIRTUAL ORBITALS FOR IRREPJR
C
       NVRTJR=VRT(IRREPJR)
C
C        DETERMINE IRREPI WHOSE DIRECT PRODUCT WITH IRREPJ GIVES IRREP
C
       IRREPI=DIRPRD(IRREPR,IRREPJR)
C
C        GET NUMBER  OF VIRTUAL ORBITALS FOR IRREPI
C
       NVRTI=VRT(IRREPI)
C
       IRREPJL=DIRPRD(IRREPI,IRREPL)
C
       NVRTJL=VRT(IRREPJL)
C
       JOFF1=1
       DO 89 IRREP=1,IRREPJL-1
        JOFF1=JOFF1+VRT(IRREP)*VRT(DIRPRD(IRREPL,IRREP))
89     CONTINUE
C
C
C        IF NVRTI OR NVRTJR OR NVRTJL EQUAL ZERO, NOTHING TO COMPUTE
C
       IF(MIN(NVRTI,NVRTJR,NVRTJL).NE.0) THEN
C
        CALL XGEMM('T','N',NVRTJL,NVRTJR,NVRTI*NUMSYR,FACT,
     &             T1(1,JOFF1),NVRTI*NUMSYL,T2(1,JOFF2),
     &             NVRTI*NUMSYR,ONE,DDVV(IOFF),NVRTJL)
       ENDIF
C
C  UPDATE OFFSETS
C
       JOFF2=JOFF2+NVRTI*NVRTJR
       IOFF=IOFF+NVRTJL*NVRTJR
90    CONTINUE
C
      RETURN
      END
