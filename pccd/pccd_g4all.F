      SUBROUTINE PCCD_G4ALL(ICORE,MAXCOR,SPCASE,IUHF)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,DISSYG,DSSYT1A,DSSYT2A,DSSYT1B,DSSYT2B,
     &        DSSYT1,DSSYT2,DSSYT3A,DSSYT3B,DSSYT4A,DSSYT4B,
     &        DISSYG2,POP,VRT  
      CHARACTER*4 SPCASE
      DIMENSION ICORE(MAXCOR)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
C
      DATA AZERO,HALF,ONE,ONEM,TWO /0.0D+0,0.5D0,1.0D+0,-1.D+0,2.D0/
C
      IOFFLISTH=0

      IF(SPCASE.EQ.'ABBA')THEN

      LISTG=118
      LISTH=156+IOFFLISTH

      DO 190 IRREP=1,NIRREP
      DISSYG=IRPDPD(IRREP,ISYTYP(1,LISTG))
      NUMSYG=IRPDPD(IRREP,ISYTYP(2,LISTG))
      I000=1
      I010=I000+IINTFP*NUMSYG*DISSYG
      IEND=I010+IINTFP*NUMSYG*DISSYG
      CALL GETLST(ICORE(I000),1,NUMSYG,1,IRREP,LISTH)
CSSS      call checksum("H4ABBA  :",ICORE(I000),DISSYG*NUMSYG)
      CALL SYMMET2(ICORE(I000),DISSYG)
      CALL PUTLST(ICORE(I000),1,NUMSYG,2,IRREP,LISTG)
 190  CONTINUE
C
      ELSEIF(SPCASE.EQ.'ABAB')THEN

       LISTG=125
       LISTH=158+IOFFLISTH

       DO 290 IRREP=1,NIRREP
       DISSYG=IRPDPD(IRREP,ISYTYP(1,LISTG))
       NUMSYG=IRPDPD(IRREP,ISYTYP(2,LISTG))
       I000=1         
       I010=I000+IINTFP*NUMSYG*DISSYG
       CALL GETLST(ICORE(I000),1,NUMSYG,1,IRREP,LISTH)
       CALL SYMMET2(ICORE(I000),DISSYG)
       CALL PUTLST(ICORE(I000),1,NUMSYG,2,IRREP,LISTG)
CSSS       call checksum("H4ABAB  :",ICORE(I000),DISSYG*NUMSYG)
290   CONTINUE

      ENDIF
      RETURN
      END
