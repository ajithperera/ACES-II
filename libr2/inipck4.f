      SUBROUTINE INIPCK4(IRREPX,IRPDPDL,IRPDPDR,LIST,IARG1,IARGX2,ISET)
C
C DRIVES THE CREATION OF MOIO POINTERS FOR SYMMETRY PACKED LISTS.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DIMENSION IRPDPDL(*),IRPDPDR(*)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      IARG2=IARGX2
      DO 10 IRREPR=1,NIRREP
       IRREPL=DIRPRD(IRREPR,IRREPX)
       CALL UPDMOI(IRPDPDR(IRREPR),IRPDPDL(IRREPL),
     &             IRREPR,LIST,IARG1,IARG2)
       IARG1=0
       IARG2=0
10    CONTINUE
C
      RETURN
      END
