      SUBROUTINE W45T3DRV(CORE,MAXCOR,IUHF,
     &                    LZOOFF,LZVOFF,LWOFF,LT3OFF,INIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C-----------------------------------------------------------------------
C     This subroutine drives the calculation of the contribution of T3
C     to the abci and mcjk Hbar elements. 
C
C     LZOOFF  offset for mcjk lists.
C     LZVOFF  offset for abci lists.
C     LWOFF   offset for ijab lists
C     LT3OFF  offset for T3 lists. T3 lists are assumed to exist and
C             the file is assumed to have been opened. LT3OFF appears
C             to be 4 for CCSDT (ie triples are 5-8), but it is 0
C             for noniterative triple excitation EOMEE methods.
C     INIT    determines whether T3 contribution is to be added to
C             existing lists (INIT = .TRUE.) or if the lists are to
C             be created and initialized. The former situation arises
C             in CCSDT energy calculations, while the latter arises in
C             CCSDR(1b), CCSDR(3), and NCCSDT-3 calculations.
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION CORE
      INTEGER MAXCOR,IUHF,LZOOFF,LZVOFF,LWOFF,LT3OFF
      LOGICAL INIT
C
      INTEGER POP,VRT,DIRPRD,SCRSIZ,DSZ,DSZEXP,DISSIZ,DISTSZ
      LOGICAL IJKEQL,IJEQL,JKEQL,NONEQL
      LOGICAL CCSDT4,CCSDT
      DIMENSION CORE(1)
      DIMENSION IADT3(8),IADW(8),IADW2(8),LEN(8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /FILES / LUOUT,MOINTS
c      COMMON /LISWI/  LWIC11,LWIC12,LWIC13,LWIC14,
c     1                LWIC15,LWIC16,LWIC17,LWIC18,
c     1                LWIC21,LWIC22,LWIC23,LWIC24,
c     1                LWIC25,LWIC26,LWIC27,LWIC28,
c     1                LWIC31,LWIC32,LWIC33,
c     1                LWIC34,LWIC35,LWIC36,
c     1                LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42
c      COMMON /GAMLIS/ LISGO1,LISGO2,LISGO3,LISGO4,LISGV1,LISGV2,
c     1                LISGV3,LISGV4
      COMMON /T3OFF / IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
      COMMON /T3IOOF/ IJKPOS(8,8,8,2),IJKLEN(36,8,4),IJKOFF(36,8,4),
     1                NCOMB(4)
      COMMON /AUXIO / DISTSZ(8,100),NDISTS(8,100),INIWRD(8,100),LNPHYR,
     1                NRECS,LUAUX
      COMMON /T3METH/ CCSDT4,CCSDT
C
      INDEX(I) = I*(I-1) / 2
C
      WRITE(LUOUT,1000)
 1000 FORMAT(' @W45T3DRV-I, Calculating triples contributions to',
     1       ' W(abci) and W(ijka). ')
C
C     -----------------------------------------------------------------
C     Block of code added for using this routine from the lambda code.
C     This is not needed when used in CCSDT calculations since the work
C     is already handled by TRPS.
C
      IRREPX = 1
C     
      CALL MKOFOO
      CALL MKOFVV
      CALL MKOFVO
c      LISGO1 = LZOOFF + 1
c      LISGO2 = LZOOFF + 2
c      LISGO3 = LZOOFF + 3
c      LISGO4 = LZOOFF + 4
c      LISGV1 = LZVOFF + 1
c      LISGV2 = LZVOFF + 2
c      LISGV3 = LZVOFF + 3
c      LISGV4 = LZVOFF + 4
C
C-----------------------------------------------------------------------
      IF(INIT)THEN
        IMODE = 0
        IF(IUHF.NE.0)THEN
          CALL INIPCK(1, 3,16,LZOOFF+1,IMODE,0,1)
          CALL INIPCK(1, 4,17,LZOOFF+2,IMODE,0,1)
          CALL INIPCK(1,14,11,LZOOFF+3,IMODE,0,1)
          CALL ZEROLIST(CORE,MAXCOR,LZOOFF+1)
          CALL ZEROLIST(CORE,MAXCOR,LZOOFF+2)
          CALL ZEROLIST(CORE,MAXCOR,LZOOFF+3)
        ENDIF
        CALL INIPCK(1,14,18,LZOOFF+4,IMODE,0,1)
        CALL ZEROLIST(CORE,MAXCOR,LZOOFF+4)
C
        IF(IUHF.NE.0)THEN
          CALL INIPCK(1, 1, 9,LZVOFF+1,IMODE,0,1)
          CALL INIPCK(1, 2,10,LZVOFF+2,IMODE,0,1)
          CALL INIPCK(1,13,18,LZVOFF+3,IMODE,0,1)
          CALL ZEROLIST(CORE,MAXCOR,LZVOFF+1)
          CALL ZEROLIST(CORE,MAXCOR,LZVOFF+2)
          CALL ZEROLIST(CORE,MAXCOR,LZVOFF+3)
        ENDIF
        CALL INIPCK(1,13,11,LZVOFF+4,IMODE,0,1)
        CALL ZEROLIST(CORE,MAXCOR,LZVOFF+4)
      ENDIF
C-----------------------------------------------------------------------
C
c      CALL INITRP
      IF(IUHF.NE.0) CALL TSPABCI(CORE,MAXCOR/IINTFP,LZVOFF+3,1)
C
C
C-----------------------------------------------------------------------
C     Reorder mcjk on disk : jk;mc --- > mc;jk.
C-----------------------------------------------------------------------
C     CALL SETOOOV(CORE(1),MAXCOR,IUHF,1)
      IF(IUHF.NE.0)THEN
        CALL TSPLST(CORE,MAXCOR,LZOOFF+1,1)
        CALL TSPLST(CORE,MAXCOR,LZOOFF+2,1)
        CALL TSPLST(CORE,MAXCOR,LZOOFF+3,1)
      ENDIF
      CALL TSPLST(CORE,MAXCOR,LZOOFF+4,1)
C-----------------------------------------------------------------------
C
      SCRSIZ = 0
      DO   5 IRREP=1,NIRREP
      MAXVV = MAX(IRPDPD(IRREP,13),IRPDPD(IRREP,19),IRPDPD(IRREP,20))
      MAXOO = MAX(IRPDPD(IRREP,14),IRPDPD(IRREP,21),IRPDPD(IRREP,22))
      MAXVAL = MAX(MAXVV,MAXOO)
      IF(SCRSIZ.LT.MAXVAL) SCRSIZ = MAXVAL
    5 CONTINUE
C
      IF(IUHF.EQ.0) GOTO 105
C
      DO 100  ISPIN=1,IUHF+1
C
      DO  90 IRPIJK=1,NIRREP
C
      DO  80   IRPK=1,NIRREP
      IF(POP(IRPK,ISPIN).EQ.0) GOTO 80
      DO  70   IRPJ=1,IRPK
      IRPJK = DIRPRD(IRPJ,IRPK)
      IRPI  = DIRPRD(IRPJK,IRPIJK)
      IRPIJ =  DIRPRD(IRPI,IRPJ)
      IRPIK =  DIRPRD(IRPI,IRPK)
C
      IF(IRPI.GT.IRPJ.OR .POP(IRPI,ISPIN).EQ.0.OR.
     1                    POP(IRPJ,ISPIN).EQ.0) GOTO 70
C
      IF(IRPI.EQ.IRPJ.AND.
     1   IRPI.EQ.IRPK.AND.POP(IRPK,ISPIN).LT.3) GOTO 70
      IF(IRPJ.EQ.IRPK.AND.POP(IRPK,ISPIN).LT.2) GOTO 70
      IF(IRPI.EQ.IRPJ.AND.POP(IRPJ,ISPIN).LT.2) GOTO 70
C
      IJKEQL = .FALSE.
      IJEQL  = .FALSE.
      JKEQL  = .FALSE.
      NONEQL = .FALSE.
      IF(IRPI.EQ.IRPJ.AND.IRPI.EQ.IRPK) IJKEQL = .TRUE.
      IF(IRPJ.EQ.IRPK.AND.IRPI.NE.IRPK)  JKEQL = .TRUE.
      IF(IRPI.EQ.IRPJ.AND.IRPI.NE.IRPK)  IJEQL = .TRUE.
      IF(IRPI.NE.IRPJ.AND.IRPI.NE.IRPK.AND.IRPJ.NE.IRPK) NONEQL = .TRUE.
C
      IF(IJKEQL)THEN
      NIJ = (POP(IRPK,ISPIN) * (POP(IRPK,ISPIN)-1))/2
      NIK = (POP(IRPK,ISPIN) * (POP(IRPK,ISPIN)-1))/2
      NJK = (POP(IRPK,ISPIN) * (POP(IRPK,ISPIN)-1))/2
      ENDIF
C
      IF(JKEQL)THEN
      NIJ = POP(IRPI,ISPIN) * POP(IRPJ,ISPIN)
      NIK = POP(IRPI,ISPIN) * POP(IRPK,ISPIN)
      NJK = (POP(IRPK,ISPIN) * (POP(IRPK,ISPIN)-1))/2
      ENDIF
C
      IF(IJEQL)THEN
      NIJ = (POP(IRPJ,ISPIN) * (POP(IRPJ,ISPIN)-1))/2
      NIK = POP(IRPI,ISPIN) * POP(IRPK,ISPIN)
      NJK = POP(IRPJ,ISPIN) * POP(IRPK,ISPIN)
      ENDIF
C
      IF(NONEQL)THEN
      NIJ = POP(IRPI,ISPIN) * POP(IRPJ,ISPIN)
      NIK = POP(IRPI,ISPIN) * POP(IRPK,ISPIN)
      NJK = POP(IRPJ,ISPIN) * POP(IRPK,ISPIN)
      ENDIF
C
C     Compute (B<C,A) lengths and addresses for given IRPIJK (=IRPABC)
C
      LENW = 0
      DO   10 IRREP =1,NIRREP
      IRPBC = DIRPRD(IRPIJK,IRREP)
      LEN(IRREP) = IRPDPD(IRPBC,ISPIN) * VRT(IRREP,ISPIN)
      LENW = LENW + LEN(IRREP)
   10 CONTINUE
C
      DO 20 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADW(IRREP) = 1
      ELSE
      IADW(IRREP) = IADW(IRREP-1) + LEN(IRREP-1)
      ENDIF
   20 CONTINUE
C
      LENT3 = DISTSZ(IRPIJK,1 + 3*(ISPIN-1))
C
      I000 = 1
      I010 = I000 + LENT3
      I020 = I010 + LENW
      NEED = IINTFP * I020
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
 1010 FORMAT(' @W45T3DRV-F, Insufficient memory. Need ',I10,/,
     1       '                                   Got  ',I10)
      STOP
      ENDIF
C
      IF(NONEQL)THEN
      KLOW  = 1
      KHIGH = POP(IRPK,ISPIN)
      JLOW  = 1
      JHIGH = POP(IRPJ,ISPIN)
      ILOW  = 1
      IHIGH = POP(IRPI,ISPIN)
      ENDIF
C
      IF(IJKEQL)THEN
      KLOW  = 3
      KHIGH = POP(IRPK,ISPIN)
      JLOW  = 2
      ILOW  = 1
      ENDIF
C
      IF(IJEQL)THEN
      KLOW  = 1
      KHIGH = POP(IRPK,ISPIN)
      JLOW  = 2
      JHIGH = POP(IRPJ,ISPIN)
      ILOW  = 1
      ENDIF
C
      IF(JKEQL)THEN
      KLOW  = 2
      KHIGH = POP(IRPK,ISPIN)
      JLOW  = 1
      ILOW  = 1
      IHIGH = POP(IRPI,ISPIN)
      ENDIF
C
      DO   60 K=KLOW,KHIGH
C
      IF(IJKEQL.OR.JKEQL) JHIGH = K-1
      DO   50 J=JLOW,JHIGH
C
      IF(IJKEQL.OR.IJEQL) IHIGH = J-1
      DO   40 I=ILOW,IHIGH
C
      CALL ZERO(CORE(I000),LENT3)
C
      IJK  = IJKPOS(IRPI,IRPJ,IRPK,1)
      IOFF = IJKOFF(IJK,IRPIJK,1+3*(ISPIN-1))
      IF(IJKEQL) IJKVAL = IOFF + ((K-1)*(K-2)*(K-3))/6 + INDEX(J-1) + I
      IF(IJEQL ) IJKVAL = IOFF + (K-1)*NIJ             + INDEX(J-1) + I
      IF(JKEQL ) IJKVAL = IOFF + (INDEX(K-1)+ J-1)*POP(IRPI,ISPIN)  + I
      IF(NONEQL) IJKVAL = IOFF + (K-1)*NIJ + (J-1)*POP(IRPI,ISPIN)  + I
C     CALL GETLIST(CORE(I000),IJKVAL,1,1,IRPIJK,1 + 4 + 3*(ISPIN-1),
C    &             IRREPX)
      CALL GETLIST(CORE(I000),IJKVAL,1,1,IRPIJK,
     &             LT3OFF + 1 + 3*(ISPIN-1),IRREPX)
C
      CALL ZERO(CORE(I010),LENW)
      CALL SYMEXPT3(CORE(I000),CORE(I010),IADW,ISPIN,IRPIJK)
C
      CALL WVT314(CORE(I010),CORE(I020),IADW,
     &             IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,
     &             ISPIN,LZVOFF,LWOFF,-1.0D+00)
      CALL WOT314(CORE(I010),CORE(I020),IADW,
     &             IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,I,J,K,ISPIN,
     &             LZOOFF,LWOFF, 1.0D+00)
C
   40 CONTINUE
   50 CONTINUE
   60 CONTINUE
C
   70 CONTINUE
   80 CONTINUE
   90 CONTINUE
  100 CONTINUE
C
  105 CONTINUE
C
      DO 210 IPASS=1,IUHF+1
C
      IF(IPASS.EQ.1)THEN
      ISPIN1 = 1
      ISPIN2 = 2
      ELSE
      ISPIN1 = 2
      ISPIN2 = 1
      ENDIF
C
      IF(IPASS.EQ.2)THEN
      CALL TSPABCI2(CORE(1),MAXCOR,LZVOFF+3,LZVOFF+4,1)
C     if general approach works, we don't need next call.
C      CALL   AT2IAB(CORE(1),0,INT2)
      ENDIF
C
c      IF(IPASS.EQ.2)THEN
      IF(IUHF.NE.0)THEN
C
C     Convert W(Ak/Ij) to W(kA/jI).
C
      I000 = 1
      I010 = I000 + SCRSIZ
      I020 = I010 + SCRSIZ
      I030 = I020 + SCRSIZ
      DO 310 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP,11)
      NDIS   = IRPDPD(IRREP,14)
      CALL GETLST(CORE(I030),1,NDIS,2,IRREP,LZOOFF+3)
      CALL SYMTR3(IRREP,VRT(1,1),POP(1,2),DISSIZ,NDIS,CORE(I030),
     1            CORE(I000),CORE(I010),CORE(I020))
      CALL SYMTR1(IRREP,POP(1,1),POP(1,2),DISSIZ,     CORE(I030),
     1            CORE(I000),CORE(I010),CORE(I020))
      CALL PUTLST(CORE(I030),1,NDIS,2,IRREP,LZOOFF+3)
  310 CONTINUE
C
      ENDIF
C
      IF(IPASS.EQ.2)THEN
C     above is really 2 if general approach works.
C
      write(6,*) ' transposing abij '
C     Convert <Ab/Ij> to <bA/jI>.
C
      I000 = 1
      I010 = I000 + SCRSIZ
      I020 = I010 + SCRSIZ
      I030 = I020 + SCRSIZ
      DO 320 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP,13)
      NDIS   = IRPDPD(IRREP,14)
      CALL GETLST(CORE(I030),1,NDIS,2,IRREP,LWOFF+3)
      CALL SYMTR3(IRREP,VRT(1,1),VRT(1,2),DISSIZ,NDIS,CORE(I030),
     1            CORE(I000),CORE(I010),CORE(I020))
      CALL SYMTR1(IRREP,POP(1,1),POP(1,2),DISSIZ,     CORE(I030),
     1            CORE(I000),CORE(I010),CORE(I020))
      CALL PUTLST(CORE(I030),1,NDIS,2,IRREP,LWOFF+3)
  320 CONTINUE
C
      ENDIF
C
      DO 200 IRPIJK=1,NIRREP
      DO 190 IRPK  =1,NIRREP
      IF(POP(IRPK,ISPIN2).EQ.0) GOTO 190
C
      KLOW  = 1
      KHIGH = POP(IRPK,ISPIN2)
C
      DO 180 IRPJ=1,NIRREP
      IF(POP(IRPJ,ISPIN1).EQ.0) GOTO 180
      IRPJK = DIRPRD(IRPJ,IRPK)
      IRPI  = DIRPRD(IRPJK,IRPIJK)
      IF(IRPI.GT.IRPJ) GOTO 180
C
      IF(POP(IRPI,ISPIN1).EQ.0) GOTO 180
C
      IJEQL  = .FALSE.
      NONEQL = .FALSE.
      IF(IRPI.EQ.IRPJ)THEN
      IJEQL  = .TRUE.
      JLOW  = 2
      JHIGH = POP(IRPJ,ISPIN1)
      ILOW  = 1
      NIJ = (POP(IRPJ,ISPIN1) * (POP(IRPJ,ISPIN1)-1))/2
      ELSE
      NONEQL = .TRUE.
      JLOW  = 1
      JHIGH = POP(IRPJ,ISPIN1)
      ILOW  = 1
      IHIGH = POP(IRPI,ISPIN1)
      NIJ =  POP(IRPI,ISPIN1) * POP(IRPJ,ISPIN1)
      ENDIF
C
      NIK = POP(IRPI,ISPIN1) * POP(IRPK,ISPIN2)
      NJK = POP(IRPJ,ISPIN1) * POP(IRPK,ISPIN2)
C
      IF(IJEQL.AND.POP(IRPJ,ISPIN1).LT.2) GOTO 180
C
      IRPIJ =  DIRPRD(IRPI,IRPJ)
      IRPIK =  DIRPRD(IRPI,IRPK)
C
      LENT3 = 0
      DO 110 IRREP=1,NIRREP
      IRPAB = DIRPRD(IRREP,IRPIJK)
      LEN(IRREP) = IRPDPD(IRPAB,ISPIN1)*VRT(IRREP,ISPIN2)
      LENT3 = LENT3 + LEN(IRREP)
  110 CONTINUE
C
      DO 120 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADT3(IRREP) = 1
      ELSE
      IADT3(IRREP) = IADT3(IRREP-1) + LEN(IRREP-1)
      ENDIF
  120 CONTINUE
C
      LENW=0
      DO 130 IRREP=1,NIRREP
      JRREP = DIRPRD(IRREP,IRPIJK)
      LEN(IRREP) = IRPDPD(JRREP,13) * VRT(IRREP,ISPIN1)
      LENW = LENW + LEN(IRREP)
  130 CONTINUE
C
      DO 140 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADW(IRREP) = 1
      ELSE
      IADW(IRREP) = IADW(IRREP-1) + LEN(IRREP-1)
      ENDIF
  140 CONTINUE
C
      I000 = 1
      I010 = I000 + LENT3
      I020 = I010 + LENW
C
C     BLOCK OF AAAA T(  ,IJ) WHERE I AND J ARE DIFFERENT SYMMETRY
C     BLOCKS. LENGTH IS IRPDPD(IRPIJ,1) * POP(IRPI,1) * POP(IRPJ,1).
C
C     LIST = 13 + ISPIN1
      LIST = LWOFF + ISPIN1
      CALL GETLST(CORE(I020),IOFFOO(IRPJ,IRPIJ,ISPIN1)+1,NIJ,1,
     1            IRPIJ,LIST)
C
C     EXPAND INTO A DIFFERENT REGION OF CORE SO WE KEEP UNTOUCHED T2
C     AS WELL AS EXPANDED ONE. LENGTH IS DSZEXP * POP(IRPI,1) * POP(IRPJ,1)
C
      DSZ   = IRPDPD(IRPIJ,   ISPIN1)
      DSZEXP= IRPDPD(IRPIJ,18+ISPIN1)
      I030 = I020 + DSZ    * NIJ
      I040 = I030 + DSZEXP * NIJ
      CALL SYMEXP2(IRPIJ,VRT(1,ISPIN1),DSZEXP,DSZ,NIJ,
     1             CORE(I030),CORE(I020))
      I050 = I040 + IRPDPD(IRPJK,13) * NJK
      I060 = I050 + IRPDPD(IRPIK,13) * NIK
C     LIST = 16
      LIST = LWOFF + 3
c      IF(IUHF.EQ.0)THEN
c      CALL GETLST(CORE(I040),IOFFOO(IRPK,IRPJK,4+ISPIN2)+1,NJK,1,
c     1            IRPJK,LIST)
c      CALL GETLST(CORE(I050),IOFFOO(IRPK,IRPIK,4+ISPIN2)+1,NIK,1,
c     1            IRPIK,LIST)
c      ELSE
      CALL GETLST(CORE(I040),IOFFOO(IRPK,IRPJK,4+ISPIN1)+1,NJK,1,
     1            IRPJK,LIST)
      CALL GETLST(CORE(I050),IOFFOO(IRPK,IRPIK,4+ISPIN1)+1,NIK,1,
     1            IRPIK,LIST)
c      ENDIF
C
      I070 = I060 + SCRSIZ
      I080 = I070 + SCRSIZ
      I090 = I080 + SCRSIZ
      I100 = I090 + LENW
      NEED = I100 * IINTFP
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
      STOP
      ENDIF
      ISTART = I090
C
      DO  170    K=KLOW,KHIGH
C
      DO  160    J=JLOW,JHIGH
C
      IF(IJEQL) IHIGH = J-1
C
      DO  150    I=ILOW,IHIGH
C
      IJK  = IJKPOS(IRPI,IRPJ,IRPK,2)
      IOFF = IJKOFF(IJK,IRPIJK,ISPIN1+1)
      IF(IJEQL ) IJKVAL = IOFF + (K-1)*NIJ + INDEX(J-1) + I
      IF(NONEQL) IJKVAL = IOFF + (K-1)*NIJ + (J-1)*POP(IRPI,ISPIN1) + I
      CALL ZERO(CORE(I000),LENT3)
C     CALL GETLIST(CORE(I000),IJKVAL,1,1,IRPIJK,ISPIN1 + 1 + 4,IRREPX)
      CALL GETLIST(CORE(I000),IJKVAL,1,1,IRPIJK,
     &             LT3OFF + ISPIN1 + 1,IRREPX)
C
      CALL ZERO(CORE(I010),LENW)
      CALL ZERO(CORE(ISTART),LENW)
      CALL SYMTRW2(CORE(I000),CORE(I010),CORE(ISTART),
     1             IADT3,IADW2,IRPIJK,ISPIN1,ISPIN2)
C
      CALL WVT323(CORE(I000),CORE(I010),CORE(ISTART),IADT3,IADW2,
     1            IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,
     1            CORE(I060),CORE(I070),CORE(I080),IUHF,
     1            CORE(I030),CORE(I040),CORE(I050),
     1            DSZEXP,IRPDPD(IRPJK,13),IRPDPD(IRPIK,13),
     1            ISPIN1,ISPIN2,LZVOFF,LWOFF,-1.0D+00)
      CALL WOT323(CORE(I000),CORE(I010),CORE(ISTART),IADT3,IADW2,
     1            IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,
     1            I,J,K,ISPIN1,ISPIN2,IUHF,LZOOFF,LWOFF, 1.0D+00)
C
  150 CONTINUE
  160 CONTINUE
  170 CONTINUE
C
  180 CONTINUE
  190 CONTINUE
  200 CONTINUE
C
      IF(IPASS.EQ.2)THEN
      CALL TSPABCI2(CORE(1),MAXCOR,LZVOFF+3,LZVOFF+4,2)
C      CALL   AT2IAB(CORE(1),1,INT2)
      ENDIF
C
C     Convert W(kA/jI) to W(Ak/Ij).
C
c      IF(IPASS.EQ.2)THEN
      IF(IUHF.NE.0)THEN
C
      I000 = 1
      I010 = I000 + SCRSIZ
      I020 = I010 + SCRSIZ
      I030 = I020 + SCRSIZ
      DO 410 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP,11)
      NDIS   = IRPDPD(IRREP,14)
      CALL GETLST(CORE(I030),1,NDIS,2,IRREP,LZOOFF+3)
      CALL SYMTR3(IRREP,POP(1,2),VRT(1,1),DISSIZ,NDIS,CORE(I030),
     1            CORE(I000),CORE(I010),CORE(I020))
      CALL SYMTR1(IRREP,POP(1,2),POP(1,1),DISSIZ,     CORE(I030),
     1            CORE(I000),CORE(I010),CORE(I020))
      CALL PUTLST(CORE(I030),1,NDIS,2,IRREP,LZOOFF+3)
  410 CONTINUE
C
      ENDIF
C
      IF(IPASS.EQ.2)THEN
C
C     Convert <bA/jI> to <Ab/Ij>.
C
      I000 = 1
      I010 = I000 + SCRSIZ
      I020 = I010 + SCRSIZ
      I030 = I020 + SCRSIZ
      DO 420 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP,13)
      NDIS   = IRPDPD(IRREP,14)
      CALL GETLST(CORE(I030),1,NDIS,2,IRREP,LWOFF+3)
      CALL SYMTR3(IRREP,VRT(1,2),VRT(1,1),DISSIZ,NDIS,CORE(I030),
     1            CORE(I000),CORE(I010),CORE(I020))
      CALL SYMTR1(IRREP,POP(1,2),POP(1,1),DISSIZ,     CORE(I030),
     1            CORE(I000),CORE(I010),CORE(I020))
      CALL PUTLST(CORE(I030),1,NDIS,2,IRREP,LWOFF+3)
  420 CONTINUE
C
      ENDIF
C
  210 CONTINUE

C
C-----------------------------------------------------------------------
C     Reorder mcjk on disk : mc;jk ---> jk;mc.
C-----------------------------------------------------------------------
C     CALL SETOOOV(CORE(1),MAXCOR,IUHF,2)
      IF(IUHF.NE.0)THEN
        CALL TSPLST(CORE,MAXCOR,LZOOFF+1,1)
        CALL TSPLST(CORE,MAXCOR,LZOOFF+2,1)
        CALL TSPLST(CORE,MAXCOR,LZOOFF+3,1)
      ENDIF
      CALL TSPLST(CORE,MAXCOR,LZOOFF+4,1)
C-----------------------------------------------------------------------
C
C     -----------------------------------------------------------------
C     Block of code added for using this routine from the lambda code.
C     This is not needed when used in CCSDT calculations since the work
C     is already handled by TRPS.
C     
      IF(IUHF.NE.0) CALL TSPABCI(CORE,MAXCOR/IINTFP,LZVOFF+3,2)
C     -----------------------------------------------------------------
C
      RETURN
      END
