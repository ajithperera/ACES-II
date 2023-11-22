      SUBROUTINE OCSYMTR1(CORE,DISSIZ,NDIS,IRREP,LIST,IOFFOV,ICODE)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION CORE(1)
      INTEGER DISSIZ,O,V,OV,VO
C
      DIMENSION IOFFOV(8,8,4)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
 
#include "icdacc.com"
 
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     1                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
      IF(ICODE.EQ.1)THEN
      WRITE(6,1010) IRREP,LIST
 1010 FORMAT(' @OCSYMTR1-I, (O,v) being transposed to (v,O).',
     1       ' Sym, list ',2I4)
      ELSE
      WRITE(6,1020) IRREP,LIST
 1020 FORMAT(' @OCSYMTR1-I, (v,O) being transposed to (O,v).',
     1       ' Sym, list ',2I4)
      ENDIF
C
      NBYTES = DISSIZ * NDIS * 8
      WRITE(6,1030) NBYTES
 1030 FORMAT(' @OCSYMTR1-I, Scratch file requires ',I12,' Bytes. ')
C
C     Open the scratch direct access file which will be used for holding
C     a symmetry block of (Ab,Ic) integrals/intermediates.
C
      LUABIC = 90
      LNPHYR = DISSIZ
c YAU - old
c      OPEN(UNIT=LUABIC,STATUS='NEW',FILE='ABIC',ACCESS='DIRECT',
c     1     FORM='UNFORMATTED',RECL=2 * idaccm * LNPHYR)
c YAU - new
      OPEN(UNIT=LUABIC,STATUS='NEW',FILE='ABIC',ACCESS='DIRECT',
     1     FORM='UNFORMATTED',RECL=LNPHYR * iintfp * idaccm)
c YAU - end
C
      I000 = 1
      I010 = I000 + DISSIZ
C
      IF(ICODE.EQ.1)THEN
C
      DO  10 OV=1,NDIS
      CALL GETLST(CORE(I000),OV,1,2,IRREP,LIST)
      CALL DUMPIT(CORE(I000),DISSIZ,LUABIC,OV)
   10 CONTINUE
C
C     The whole list has been read and has been put on ABIC in original
C     order. Now read it and put it back on MOINTS in transposed order.
C
      IRPOV = IRREP
      IRPVO = IRREP
      DO  50 IRPV=1,NIRREP
      IRPO = DIRPRD(IRPV,IRPOV)
C
      IF(VRT(IRPV,2).EQ.0.OR.POP(IRPO,1).EQ.0) GOTO 50
C
      DO  40    V=1,VRT(IRPV,2)
      DO  30    O=1,POP(IRPO,1)
C
      OV = IOFFOV(IRPV,IRPOV,3) + (V-1) * POP(IRPO,1) + O
      VO = IOFFVO(IRPO,IRPVO,3) + (O-1) * VRT(IRPV,2) + V
C
      CALL PULLIT(CORE(I000),DISSIZ,LUABIC,OV)
      CALL PUTLST(CORE(I000),VO,1,2,IRPVO,LIST)
   30 CONTINUE
   40 CONTINUE
   50 CONTINUE
C
      ENDIF
C
      IF(ICODE.EQ.2)THEN
C
      DO 110 VO=1,NDIS
      CALL GETLST(CORE(I000),VO,1,2,IRREP,LIST)
      CALL DUMPIT(CORE(I000),DISSIZ,LUABIC,VO)
  110 CONTINUE
C
C     The whole list has been read and has been put on ABIC in original
C     order. Now read it and put it back on MOINTS in transposed order.
C
      IRPOV = IRREP
      IRPVO = IRREP
      DO 150 IRPO=1,NIRREP
      IRPV = DIRPRD(IRPO,IRPVO)
C
      IF(VRT(IRPV,2).EQ.0.OR.POP(IRPO,1).EQ.0) GOTO 150
C
      DO 140    O=1,POP(IRPO,1)
      DO 130    V=1,VRT(IRPV,2)
C
      OV = IOFFOV(IRPV,IRPOV,3) + (V-1) * POP(IRPO,1) + O
      VO = IOFFVO(IRPO,IRPVO,3) + (O-1) * VRT(IRPV,2) + V
C
      CALL PULLIT(CORE(I000),DISSIZ,LUABIC,VO)
      CALL PUTLST(CORE(I000),OV,1,2,IRPOV,LIST)
  130 CONTINUE
  140 CONTINUE
  150 CONTINUE
C
      ENDIF
C
C     Close and delete the scratch file.
C
      CLOSE(UNIT=LUABIC,STATUS='DELETE')
      RETURN
      END
