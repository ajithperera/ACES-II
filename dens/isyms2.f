c------------------------------------------------------------------
      INTEGER FUNCTION isyms2(ITYPL,ITYPR)
C
      IMPLICIT INTEGER (A-Z)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                 DIRPRD(8,8)
      COMMON /SYMPOP2/ IRPDPD(8,22)
      isyms2=0
      DO 10 IRREP=1,NIRREP
       isyms2=isyms2+IRPDPD(IRREP,ITYPL)*IRPDPD(IRREP,ITYPR)
10    CONTINUE
      RETURN
      END
