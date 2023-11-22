
      SUBROUTINE SYMT2AB(IRREPX,SCR,TMP,MAXCOR)
C
C THIS ROUTINE ACCEPTS THE T1(AA) AND T2(ABAB) VECTORS AND
C RETURNS THE SPIN-ADAPTED T VECTOR [2*T1(AA) ; 2 T(IJ,AB) - T(IJ,BA)].
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT,DISSIZ
      DIMENSION SCR(MAXCOR),TMP(*)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
C
      LENC1A=IRPDPD(IRREPX,9)
      LENC2AB=IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
      I0C1A=1
      I0C2AB=I0C1A+IRPDPD(IRREPX,9)
      ITOP =I0C2AB+LENC2AB
C
      HALF=0.5D0
C
      IOFF=I0C2AB
      ITMP1=1
      DO 10 IRREPR=1,NIRREP
       IRREPL=DIRPRD(IRREPR,IRREPX)
       DISSIZ=IRPDPD(IRREPL,ISYTYP(1,16))
       NUMDIS=IRPDPD(IRREPR,ISYTYP(2,16))
       ITMP2=ITMP1+MAX(DISSIZ,NUMDIS)
       ITMP3=ITMP2+MAX(DISSIZ,NUMDIS)
       CALL SYMRHF3(IRREPL,IRREPR,VRT(1,1),POP(1,1),DISSIZ,
     &              SCR(IOFF),TMP(ITMP1),TMP(ITMP2),TMP(ITMP3))
       CALL SSCAL(DISSIZ*NUMDIS,HALF,SCR(IOFF),1)
       IOFF=IOFF+DISSIZ*NUMDIS
10    CONTINUE
C
      RETURN
      END