      

 
      SUBROUTINE CHECKGAM(ICORE,LISTW,LISTG,FACT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,DISSYW
      DIMENSION ICORE(*)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON/ADD/SUM
CSSS      RETURN
      E=0.0D+0
      E1=0.0D+0
      E2=0.0D+0
      sum = 0.0D0
      DO 1000 IRREP=1,NIRREP
      NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
      DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
      IOFFW=1
      IOFFW2=1+NUMSYW*DISSYW*IINTFP
      CALL GETLST(ICORE(IOFFW),1,NUMSYW,1,IRREP,LISTW)
      CALL GETLST(ICORE(IOFFW2),1,NUMSYW,2,IRREP,LISTG)
CSSS      call output(icore(ioffw2),1,dissyw,1,numsyw,dissyw,numsyw,1)
      E=E+SDOT(NUMSYW*DISSYW,ICORE(IOFFW),1,ICORE(IOFFW2),1)
      E1=E1+SDOT(NUMSYW*DISSYW,ICORE(IOFFW),1,ICORE(IOFFW),1)
      E2=E2+SDOT(NUMSYW*DISSYW,ICORE(IOFFW2),1,ICORE(IOFFW2),1)
1000  CONTINUE
      sum=sum+FACT*e
      write(6,500)Fact*e,e1,e2
      write(6,501)sum
      return
500   format(' Energy contribution : ',3(F13.10,1X))
501   format(' Cumulative energy   : ',F20.10)
      end
