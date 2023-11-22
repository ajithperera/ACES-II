      SUBROUTINE ACTDNMBL(CORE,idimcore,MAXCOR,ISPIN,STRING,FACTOR)
C
C     THIS ROUTINE PALACES FACTOR INTO THE DENOMINATOR DESCRIBED
C     BY STRING IN ORDER TO EXCLUDE CERTAIN T2 AMPLITUDES
C     FROM THE CALCULATIONS
C     
CEND
C
C CODED PS SEPT/92
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION CORE,FACTOR,ZILCH
      DIMENSION CORE(idimcore)
C
      CHARACTER*4 STRING
      CHARACTER*2 TYPLEFT,TYPRGHT
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/ IFLAGS(100)
C
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
C
      COMMON /FSSYMPOP/ FSDPDAN(8,22),FSDPDNA(8,22),FSDPDAA(8,22),
     &                  FSDPDIN(8,22),FSDPDNI(8,22),FSDPDII(8,22),
     &                  FSDPDAI(8,22),FSDPDIA(8,22)
C
      COMMON /FSSYM/ POPA(8,2),POPI(8,2),VRTA(8,2),VRTI(8,2),
     &               NTAN(2),NTNA(2),NTAA(2),
     &               NTIN(2),NTNI(2),NTII(2),
     &               NTAI(2),NTIA(2),
     &               NF1AN(2),NF1AA(2),NF1IN(2),NF1II(2),NF1AI(2),
     &               NF2AN(2),NF2AA(2),NF2IN(2),NF2II(2),NF2AI(2)
C
      DATA ZILCH/0.D0/
      DATA LISTD0/47/
C
      IF(ISPIN.NE.3) THEN
         WRITE(6,1100) ISPIN
 1100    FORMAT(T3,'@ACTDNMBL-E ISPIN=',A,
     $      ' Sorry, only ISPIN=3 is allowed')
         CALL ERREX
      ENDIF
C
      LISTD=LISTD0+ISPIN
      TYPLEFT=STRING(1:2)
      TYPRGHT=STRING(3:4)
C
      DO 1 IRREP=1,NIRREP
C
C   DO FOR PARTICLE LABELES (LEFT)
C
      IF(TYPLEFT.EQ.'NA')THEN
         DISSYT=FSDPDNA(IRREP,ISYTYP(1,LISTD))
      ELSE IF(TYPLEFT.EQ.'NI')THEN
         DISSYT=FSDPDNI(IRREP,ISYTYP(1,LISTD))
      ELSE IF(TYPLEFT.EQ.'AN')THEN
         DISSYT=FSDPDAN(IRREP,ISYTYP(1,LISTD))
      ELSE IF(TYPLEFT.EQ.'IN')THEN
         DISSYT=FSDPDIN(IRREP,ISYTYP(1,LISTD))
      ELSE IF(TYPLEFT.EQ.'AA')THEN
         DISSYT=FSDPDAA(IRREP,ISYTYP(1,LISTD))
      ELSE IF(TYPLEFT.EQ.'AI')THEN
         DISSYT=FSDPDAI(IRREP,ISYTYP(1,LISTD))
      ELSE IF(TYPLEFT.EQ.'IA')THEN
         DISSYT=FSDPDIA(IRREP,ISYTYP(1,LISTD))
      ELSE IF(TYPLEFT.EQ.'II')THEN
         DISSYT=FSDPDII(IRREP,ISYTYP(1,LISTD))
      ELSE IF(TYPLEFT.EQ.'NN')THEN
         DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTD))
      ELSE
         WRITE(6,1200) TYPLEFT
 1200    FORMAT(T3,'@ACTDNMBL-E TYPLEFT=',A,' does not make sense')
         CALL ERREX
      ENDIF
C
C   DO FOR HOLE LABELES (RIGHT)
C
      IF(TYPRGHT.EQ.'NA')THEN
         NUMSYT=FSDPDNA(IRREP,ISYTYP(2,LISTD))
      ELSE IF(TYPRGHT.EQ.'NI')THEN
         NUMSYT=FSDPDNI(IRREP,ISYTYP(2,LISTD))
      ELSE IF(TYPRGHT.EQ.'AN')THEN
         NUMSYT=FSDPDAN(IRREP,ISYTYP(2,LISTD))
      ELSE IF(TYPRGHT.EQ.'IN')THEN
         NUMSYT=FSDPDIN(IRREP,ISYTYP(2,LISTD))
      ELSE IF(TYPRGHT.EQ.'AA')THEN
         NUMSYT=FSDPDAA(IRREP,ISYTYP(2,LISTD))
      ELSE IF(TYPRGHT.EQ.'AI')THEN
         NUMSYT=FSDPDAI(IRREP,ISYTYP(2,LISTD))
      ELSE IF(TYPRGHT.EQ.'IA')THEN
         NUMSYT=FSDPDIA(IRREP,ISYTYP(2,LISTD))
      ELSE IF(TYPRGHT.EQ.'II')THEN
         NUMSYT=FSDPDII(IRREP,ISYTYP(2,LISTD))
      ELSE IF(TYPRGHT.EQ.'NN')THEN
         NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTD))
      ELSE
         WRITE(6,1300) TYPRGHT
 1300    FORMAT(T3,'@ACTDNMBL-E TYPRGHT=',A,' does not make sense')
         CALL ERREX
      ENDIF
C
      I000=1
      I010=I000+NUMSYT*DISSYT*IINTFP
      IF(NUMSYT*DISSYT.EQ.0) GOTO 1
      IF(I010.GT.MAXCOR) CALL INSMEM('ACTDNBL',I010,MAXCOR)
C
      DO 10 I=1,NUMSYT*DISSYT
         CORE(I)=FACTOR
 10   CONTINUE
C      
      CALL FSPUT(CORE(I000),1,NUMSYT,1,IRREP,LISTD+16,STRING)
C
      IF(FACTOR.NE.ZILCH) THEN
         DO 11 I=1,NUMSYT*DISSYT
            CORE(I)=1.D0/FACTOR
 11      CONTINUE
C      
         CALL FSPUT(CORE(I000),1,NUMSYT,1,IRREP,LISTD,STRING)
      ELSE
         WRITE(6,1400) FACTOR,LISTD+16
 1400    FORMAT(T3,'@ACTDNMBL-W FACTOR=',F10.6,' is zero, list=',i3,
     $      'is not modified')
      ENDIF
 1    CONTINUE
      RETURN
      END
