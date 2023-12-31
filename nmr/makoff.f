      SUBROUTINE MAKOFF(IOFFU,IOFFS1,IOFFS2,IRREPX)
C
C  THIS ROUTINE COMPUTES OFFSETS FOR THE U AND S MATRICES
C
CEND
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER POP,VRT,DIRPRD
C
      DIMENSION IOFFU(8,2),IOFFS1(8,2),IOFFS2(8,2) 
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
C
      IOFFU(1,1)=1
      IOFFU(1,2)=1
      IOFFS1(1,1)=1
      IOFFS1(1,2)=1
      IOFFS2(1,1)=1
      IOFFS2(1,2)=1
      DO 10 IRREPR=1,NIRREP-1
       IRREPL=DIRPRD(IRREPX,IRREPR)
       IOFFU(IRREPR+1,1)=IOFFU(IRREPR,1)+POP(IRREPR,1)*VRT(IRREPL,1)
       IOFFU(IRREPR+1,2)=IOFFU(IRREPR,2)+POP(IRREPR,2)*VRT(IRREPL,2)
       IOFFS1(IRREPR+1,1)=IOFFS1(IRREPR,1)+POP(IRREPR,1)*POP(IRREPL,1)
       IOFFS1(IRREPR+1,2)=IOFFS1(IRREPR,2)+POP(IRREPR,2)*POP(IRREPL,2)
       IOFFS2(IRREPR+1,1)=IOFFS2(IRREPR,1)+VRT(IRREPR,1)*VRT(IRREPL,1)
       IOFFS2(IRREPR+1,2)=IOFFS2(IRREPR,2)+VRT(IRREPR,2)*VRT(IRREPL,2)
10    CONTINUE
C
      RETURN
      END
