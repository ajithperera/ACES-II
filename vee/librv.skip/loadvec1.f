      SUBROUTINE LOADVEC1(IRREPX,SCR,MAXCOR,IUHF,LIST1,IOF1,
     &                    LIST2,LENTOT,SPINAD)
C
C THIS ROUTINE LOADS A VECTOR INTO CORE, ORDERED AS FOLLOWS 
C
C   SINGLE(AA)-SINGLE(BB)-DOUBLE(AA)-DOUBLE(BB)-DOUBLE(AB)  [UHF]
C
C   SINGLE(AA)-DOUBLE(AB)                                   [RHF]
C
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT,DISSIZ
      LOGICAL SPINAD
      DIMENSION SCR(MAXCOR)
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
C
      LENC1A=IRPDPD(IRREPX,9)
      LENC2AB=IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
      I0C1A=1
      CALL GETLST(SCR(I0C1A),1,1,1,1+IOF1,LIST1)
      IF(IUHF.NE.0)THEN
       LENC1B=IRPDPD(IRREPX,10)
       LENC2AA=IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
       LENC2BB=IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
       I0C1B=I0C1A+LENC1A
       I0C2AB=I0C1B+LENC1B
       I0C2BB=I0C2AB+LENC2AB
       I0C2AA=I0C2BB+LENC2BB 
       ITOP  =I0C2AA+LENC2AA
       CALL GETLST(SCR(I0C1B),1,1,1,2+IOF1,LIST1)
       CALL GETALL(SCR(I0C2AA),1,IRREPX,LIST2+1)
       CALL GETALL(SCR(I0C2BB),1,IRREPX,LIST2+2)
      ELSE
       I0C2AB=I0C1A+LENC1A
       ITOP  =I0C2AB+LENC2AB
      ENDIF
      CALL GETALL (SCR(I0C2AB),1,IRREPX,LIST2+3)
C
      IF(SPINAD)THEN
       IOFF=I0C2AB
       ITMP1=ITOP
       DO 10 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPX)
        DISSIZ=IRPDPD(IRREPL,ISYTYP(1,16))
        NUMDIS=IRPDPD(IRREPR,ISYTYP(2,16))
        ITMP2=ITMP1+MAX(DISSIZ,NUMDIS)
        CALL SPINAD3(IRREPL,VRT(1,1),DISSIZ,NUMDIS,SCR(IOFF),SCR(ITMP1),
     &                SCR(ITMP2))
        IOFF=IOFF+DISSIZ*NUMDIS
10     CONTINUE
       CALL SSCAL(LENC1A,2.0D0,SCR(I0C1A),1)
      ENDIF
C
      LENTOT=ITOP-1
C
      RETURN
      END
