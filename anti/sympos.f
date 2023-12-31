      SUBROUTINE SYMPOS(TYPE,POP1,POP2,IRREP,IRDPD,OFFSET)
      IMPLICIT INTEGER (A-Z)
      DIMENSION POP1(8),POP2(8)
      CHARACTER*1 TYPE
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
C
      NNM1O2(I)=(I*(I-1))/2
C
      OFFSET=0
      IF(TYPE.EQ.'F')THEN
       DO 10 IRREPX=1,IRREP-1
        IRREPY=DIRPRD(IRREPX,IRDPD)
        OFFSET=OFFSET+POP1(IRREPX)*POP2(IRREPY)
10     CONTINUE
      ELSEIF(TYPE.EQ.'T')THEN
       DO 20 IRREPX=1,IRREP-1
        IRREPY=DIRPRD(IRREPX,IRDPD)
        IF(IRREPX.LT.IRREPY)THEN
         OFFSET=OFFSET+POP1(IRREPX)*POP2(IRREPY)
        ELSEIF(IRREPY.LT.IRREPX)THEN
         OFFSET=OFFSET+POP1(IRREPY)*POP2(IRREPX)
        ELSEIF(IRREPY.EQ.IRREPX)THEN
         OFFSET=OFFSET+NNM1O2(POP1(IRREPX))
        ENDIF
20     CONTINUE
      ENDIF
      RETURN
      END
