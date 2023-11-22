      SUBROUTINE T3WT323R(D3T3,D3T3EXP,CORE,MAXCOR,ISPIN1,ISPIN2,
     1                    IADDT,IADDTX,IRPI,IRPJ,IRPK,I,J,K,IRPIJK,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER VRT,POP,DIRPRD,DISSIZ,DISTSZ
      LOGICAL NONEQL,
     1        IJMEQL,IMEQL,JMEQL,IJEQL,
     1               MIEQL,MJEQL
      DIMENSION D3T3(1),D3T3EXP(1),CORE(1)
      DIMENSION IADDT(8),IADDTX(8)
      DIMENSION IADT(8),LENT(8),IADTX(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NT1(2),NFMI(2),NFEA(2)
      COMMON /FLAGS/  IFLAGS(100)
      EQUIVALENCE(ICLLVL,IFLAGS( 2))
      EQUIVALENCE(IDRLVL,IFLAGS( 3))
      COMMON /FILES/  LUOUT,MOINTS
C
      COMMON /AUXIO/  DISTSZ(8,100),NDISTS(8,100),INIWRD(8,100),LNPHYR,
     1                NRECS,LUAUX
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
      COMMON /T3IOOF/ IJKPOS(8,8,8,2),IJKLEN(36,8,4),IJKOFF(36,8,4),
     1                NCOMB(4)
C
      INDEX(I) = I*(I-1)/2
C
      DO  100 IRPM=1,NIRREP
C
      IF(POP(IRPM,ISPIN1).EQ.0) GOTO 100
      IF(POP(IRPM,ISPIN1).LT.2.AND.IRPM.EQ.IRPJ) GOTO 100
C
      IF(IRPM.LE.IRPJ) MJKPOS = IJKPOS(IRPM,IRPJ,IRPK,2)
      IF(IRPM.GT.IRPJ) MJKPOS = IJKPOS(IRPJ,IRPM,IRPK,2)
C
      MJEQL  = .FALSE.
      MJEQL  = IRPM.EQ.IRPJ
C
      IRPJK  = DIRPRD(IRPJ,IRPK)
      IRPMJK = DIRPRD(IRPM,IRPJK)
C
      IOFF   = IJKOFF(MJKPOS,IRPMJK,ISPIN1 + 1)
C
      LN    = DISTSZ(IRPMJK,ISPIN1 + 1 + 4)
      LNEXP = 0
      DO 10 IRREP=1,NIRREP
      JRREP = DIRPRD(IRREP,IRPMJK)
      LNEXP = LNEXP + IRPDPD(JRREP,13) * VRT(IRREP,ISPIN1)
      LENT(IRREP) = IRPDPD(JRREP,ISPIN1) * VRT(IRREP,ISPIN2)
   10 CONTINUE
C
      DO 20 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADT(IRREP) = 1
      ELSE
      IADT(IRREP) = IADT(IRREP-1) + LENT(IRREP-1)
      ENDIF
   20 CONTINUE
C
C     I000 -  W(C,E)
C     I010 - T3(A<B,e)
C     I020 - T3(B,e,A) (labelled by A)
C
      IRPMI = DIRPRD(IRPI,IRPM)
      I000 = 1
      I010 = I000 + MAX(IRPDPD(IRPMI,18 + ISPIN1),
     1                  IRPDPD(IRPMI,18 + ISPIN2))
      I020 = I010 + LN
      I030 = I020 + LNEXP
      I040 = I030 + LNEXP
      NEED = IINTFP * I040
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
 1010 FORMAT(' @T3WT323R-F, Insufficient memory. Need ',I10,' Got ',I10)
      CALL INSMEM('T3WT323R',NEED,MAXCOR)
      ENDIF
C
      DO   90    M=1,POP(IRPM,ISPIN1)
C
      SIGN =  1.0D+00
C
      IF(MJEQL.AND.M.EQ.J) GOTO 90
C
      IF(MJEQL)THEN
       NMJ = (POP(IRPM,ISPIN1) * (POP(IRPM,ISPIN1)-1))/2
       IF(M.LT.J)THEN
       MJ = INDEX(J-1) + M
       ELSE
       MJ = INDEX(M-1) + J
       SIGN = - SIGN
       ENDIF
      ELSE
       NMJ = POP(IRPM,ISPIN1) * POP(IRPJ,ISPIN1)
       IF(IRPJ.GT.IRPM)THEN
       MJ = (J-1)*POP(IRPM,ISPIN1) + M
       ELSE
       MJ = (M-1)*POP(IRPJ,ISPIN1) + J
       SIGN = - SIGN
       ENDIF
      ENDIF
      MJKVAL = IOFF + (K-1)*NMJ + MJ
C
      CALL GETLIST(CORE(I010),MJKVAL,1,1,IRPMJK,ISPIN1 + 1 + 4)
C
      MI = IOFFOO(IRPI,IRPMI,2+ISPIN1) + (I-1)*POP(IRPM,ISPIN1) + M
C
C     Direct summation into D3T3 : D3T3(A<B,c) = T3(A<B,e) * W(e,c)
C                                  D3T3(a<b,C) = T3(a<b,E) * W(E,C)
C
      IF(IUHF.EQ.0)THEN
      LISTW = 58
      ELSE
      LISTW = 60 - ISPIN1
      ENDIF
      CALL GETLST(CORE(I000),MI,1,2,IRPMI,LISTW)
C
      IOFFT = I010
      DO 50 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN2).EQ.0) GOTO 50
C
      IRPC  = DIRPRD(IRPE,IRPMI)
      IRPCE = IRPMI
      IOFFW = IOFFVV(IRPC,IRPCE,ISPIN2 + 2) + 1
      IRPAB = DIRPRD(IRPE,IRPMJK)
      DISSIZ = IRPDPD(IRPAB,ISPIN1)
      NDIS   = VRT(IRPC,ISPIN2)
      NSUM   = VRT(IRPE,ISPIN2)
      IOFFDT = IADDT(IRPC)
C
      CALL XGEMM('N','N',DISSIZ,NDIS,NSUM,SIGN,
     1           CORE(IOFFT),DISSIZ,CORE(IOFFW),NSUM,1.0D+00,
     1           D3T3(IOFFDT),DISSIZ)
C
      IOFFT = IOFFT + DISSIZ * NSUM
   50 CONTINUE
C
C     Summation into expanded D3T3 : D3T3EXP(A,Bc) = W(A,E) * T3EXP(Bc,E)
C                                    D3T3EXP(a,bC) = W(a,e) * T3EXP(bC,e)
C
      LISTW = 53 + ISPIN1
      CALL GETLST(CORE(I000),MI,1,2,IRPMI,LISTW)
C
      CALL SYMTRW2(CORE(I010),CORE(I020),CORE(I030),IADT,IADTX,IRPMJK,
     1             ISPIN1,ISPIN2)
C
      IOFFT = I020
      DO 60 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN1).EQ.0) GOTO 60
C
      IRPA  = DIRPRD(IRPE,IRPMI)
      IRPAE = IRPMI
      IOFFW = IOFFVV(IRPA,IRPAE,ISPIN1 + 2) + 1
      IRPBC = DIRPRD(IRPE,IRPMJK)
      DISSIZ = VRT(IRPA,ISPIN1)
      NDIS   = IRPDPD(IRPBC,13)
      NSUM   = VRT(IRPE,ISPIN1)
      IOFFDT = IADDTX(IRPBC)
C
      CALL XGEMM('T','T',DISSIZ,NDIS,NSUM,SIGN,
     1           CORE(IOFFW),NSUM,CORE(IOFFT),NDIS,1.0D+00,
     1           D3T3EXP(IOFFDT),DISSIZ)
C
      IOFFT = IOFFT + NDIS * NSUM
   60 CONTINUE
C
   90 CONTINUE
  100 CONTINUE
C
      DO  200 IRPM=1,NIRREP
C
      IF(POP(IRPM,ISPIN1).EQ.0) GOTO 200
      IF(POP(IRPM,ISPIN1).LT.2.AND.IRPM.EQ.IRPI) GOTO 200
C
      IF(IRPM.LE.IRPI) MIKPOS = IJKPOS(IRPM,IRPI,IRPK,2)
      IF(IRPM.GT.IRPI) MIKPOS = IJKPOS(IRPI,IRPM,IRPK,2)
C
      MIEQL  = .FALSE.
      MIEQL  = IRPM.EQ.IRPI
C
      IRPIK  = DIRPRD(IRPI,IRPK)
      IRPMIK = DIRPRD(IRPM,IRPIK)
C
      IOFF   = IJKOFF(MIKPOS,IRPMIK,ISPIN1 + 1)
C
      LN    = DISTSZ(IRPMIK,ISPIN1 + 1 + 4)
      LNEXP = 0
      DO 110 IRREP=1,NIRREP
      JRREP = DIRPRD(IRREP,IRPMIK)
      LNEXP = LNEXP + IRPDPD(JRREP,13) * VRT(IRREP,ISPIN1)
      LENT(IRREP) = IRPDPD(JRREP,ISPIN1) * VRT(IRREP,ISPIN2)
  110 CONTINUE
C
      DO 120 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADT(IRREP) = 1
      ELSE
      IADT(IRREP) = IADT(IRREP-1) + LENT(IRREP-1)
      ENDIF
  120 CONTINUE
C
C     I000 -  W(C,E)
C     I010 - T3(A<B,e)
C     I020 - T3(B,e,A) (labelled by A)
C
      IRPMJ = DIRPRD(IRPJ,IRPM)
      I000 = 1
      I010 = I000 + MAX(IRPDPD(IRPMJ,18 + ISPIN1),
     1                  IRPDPD(IRPMJ,18 + ISPIN2))
      I020 = I010 + LN
      I030 = I020 + LNEXP
      I040 = I030 + LNEXP
      NEED = IINTFP * I040
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
      CALL INSMEM('T3WT323R',NEED,MAXCOR)
      ENDIF
C
      DO  190    M=1,POP(IRPM,ISPIN1)
C
      SIGN = -1.0D+00
C
      IF(MIEQL.AND.M.EQ.I) GOTO 190
C
      IF(MIEQL)THEN
       NMI = (POP(IRPM,ISPIN1) * (POP(IRPM,ISPIN1)-1))/2
       IF(M.LT.I)THEN
       MI = INDEX(I-1) + M
       ELSE
       MI = INDEX(M-1) + I
       SIGN = - SIGN
       ENDIF
      ELSE
       NMI = POP(IRPM,ISPIN1) * POP(IRPI,ISPIN1)
       IF(IRPI.GT.IRPM)THEN
       MI = (I-1)*POP(IRPM,ISPIN1) + M
       ELSE
       MI = (M-1)*POP(IRPI,ISPIN1) + I
       SIGN = - SIGN
       ENDIF
      ENDIF
      MIKVAL = IOFF + (K-1)*NMI + MI
C
      CALL GETLIST(CORE(I010),MIKVAL,1,1,IRPMIK,ISPIN1 + 1 + 4)
C
      MJ = IOFFOO(IRPJ,IRPMJ,2+ISPIN1) + (J-1)*POP(IRPM,ISPIN1) + M
C
C     Direct summation into D3T3 : D3T3(A<B,c) = T3(A<B,e) * W(e,c)
C                                  D3T3(a<b,C) = T3(a<b,E) * W(E,C)
C
      IF(IUHF.EQ.0)THEN
      LISTW = 58
      ELSE
      LISTW = 60 - ISPIN1
      ENDIF
      CALL GETLST(CORE(I000),MJ,1,2,IRPMJ,LISTW)
C
      IOFFT = I010
      DO 150 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN2).EQ.0) GOTO 150
C
      IRPC  = DIRPRD(IRPE,IRPMJ)
      IRPCE = IRPMJ
      IOFFW = IOFFVV(IRPC,IRPCE,ISPIN2 + 2) + 1
      IRPAB = DIRPRD(IRPE,IRPMIK)
      DISSIZ = IRPDPD(IRPAB,ISPIN1)
      NDIS   = VRT(IRPC,ISPIN2)
      NSUM   = VRT(IRPE,ISPIN2)
      IOFFDT = IADDT(IRPC)
C
      CALL XGEMM('N','N',DISSIZ,NDIS,NSUM,SIGN,
     1           CORE(IOFFT),DISSIZ,CORE(IOFFW),NSUM,1.0D+00,
     1           D3T3(IOFFDT),DISSIZ)
C
      IOFFT = IOFFT + DISSIZ * NSUM
  150 CONTINUE
C
C     Summation into expanded D3T3 : D3T3EXP(A,Bc) = W(A,E) * T3EXP(Bc,E)
C                                    D3T3EXP(a,bC) = W(a,e) * T3EXP(bC,e)
C
      LISTW = 53 + ISPIN1
      CALL GETLST(CORE(I000),MJ,1,2,IRPMJ,LISTW)
C
      CALL SYMTRW2(CORE(I010),CORE(I020),CORE(I030),IADT,IADTX,IRPMIK,
     1             ISPIN1,ISPIN2)
C
      IOFFT = I020
      DO 160 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN1).EQ.0) GOTO 160
C
      IRPA  = DIRPRD(IRPE,IRPMJ)
      IRPAE = IRPMJ
      IOFFW = IOFFVV(IRPA,IRPAE,ISPIN1 + 2) + 1
      IRPBC = DIRPRD(IRPE,IRPMIK)
      DISSIZ = VRT(IRPA,ISPIN1)
      NDIS   = IRPDPD(IRPBC,13)
      NSUM   = VRT(IRPE,ISPIN1)
      IOFFDT = IADDTX(IRPBC)
C
      CALL XGEMM('T','T',DISSIZ,NDIS,NSUM,SIGN,
     1           CORE(IOFFW),NSUM,CORE(IOFFT),NDIS,1.0D+00,
     1           D3T3EXP(IOFFDT),DISSIZ)
C
      IOFFT = IOFFT + NDIS * NSUM
  160 CONTINUE
C
  190 CONTINUE
  200 CONTINUE
C
      SIGN = 1.0D+00
      DO  300 IRPM=1,NIRREP
C
      IF(POP(IRPM,ISPIN2).EQ.0) GOTO 300
C
      IJMPOS = IJKPOS(IRPI,IRPJ,IRPM,2)
C
      IRPIJ  = DIRPRD(IRPI,IRPJ)
      IRPIJM = DIRPRD(IRPM,IRPIJ)
C
      IOFF   = IJKOFF(IJMPOS,IRPIJM,ISPIN1 + 1)
C
      LN    = DISTSZ(IRPIJM,ISPIN1 + 1 + 4)
      LNEXP = 0
      DO 210 IRREP=1,NIRREP
      JRREP = DIRPRD(IRREP,IRPIJM)
      LNEXP = LNEXP + IRPDPD(JRREP,13) * VRT(IRREP,ISPIN1)
      LENT(IRREP) = IRPDPD(JRREP,ISPIN1) * VRT(IRREP,ISPIN2)
  210 CONTINUE
C
      DO 220 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADT(IRREP) = 1
      ELSE
      IADT(IRREP) = IADT(IRREP-1) + LENT(IRREP-1)
      ENDIF
  220 CONTINUE
C
C     I000 -  W(C,E)
C     I010 - T3(A<B,e)
C     I020 - T3(B,e,A) (labelled by A)
C
      IRPMK = DIRPRD(IRPK,IRPM)
      I000 = 1
      I010 = I000 + MAX(IRPDPD(IRPMK,18 + ISPIN1),
     1                  IRPDPD(IRPMK,18 + ISPIN2))
      I020 = I010 + LN
      I030 = I020 + LNEXP
      I040 = I030 + LNEXP
      NEED = IINTFP * I040
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
      CALL INSMEM('T3WT323R',NEED,MAXCOR)
      ENDIF
C
      DO  290    M=1,POP(IRPM,ISPIN2)
C
      IF(IRPI.EQ.IRPJ)THEN
       NIJ = (POP(IRPJ,ISPIN1) * (POP(IRPJ,ISPIN1)-1))/2
       IJ  = INDEX(J-1) + I
      ELSE
       NIJ =  POP(IRPI,ISPIN1) * POP(IRPJ,ISPIN1)
       IJ  =  (J-1)*POP(IRPI,ISPIN1) + I
      ENDIF
      IJMVAL = IOFF + (M-1)*NIJ + IJ
C
      CALL GETLIST(CORE(I010),IJMVAL,1,1,IRPIJM,ISPIN1 + 1 + 4)
C
      MK = IOFFOO(IRPK,IRPMK,2+ISPIN2) + (K-1)*POP(IRPM,ISPIN2) + M
C
C     Direct summation into D3T3 : D3T3(A<B,c) = T3(A<B,e) * W(e,c)
C                                  D3T3(a<b,C) = T3(a<b,E) * W(E,C)
C
      IF(IUHF.EQ.0)THEN
      LISTW = 54
      ELSE
      LISTW = 56 - ISPIN1
      ENDIF
      CALL GETLST(CORE(I000),MK,1,2,IRPMK,LISTW)
C
      IOFFT = I010
      DO 250 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN2).EQ.0) GOTO 250
C
      IRPC  = DIRPRD(IRPE,IRPMK)
      IRPCE = IRPMK
      IOFFW = IOFFVV(IRPC,IRPCE,ISPIN2 + 2) + 1
      IRPAB = DIRPRD(IRPE,IRPIJM)
      DISSIZ = IRPDPD(IRPAB,ISPIN1)
      NDIS   = VRT(IRPC,ISPIN2)
      NSUM   = VRT(IRPE,ISPIN2)
      IOFFDT = IADDT(IRPC)
C
      CALL XGEMM('N','N',DISSIZ,NDIS,NSUM,SIGN,
     1           CORE(IOFFT),DISSIZ,CORE(IOFFW),NSUM,1.0D+00,
     1           D3T3(IOFFDT),DISSIZ)
C
      IOFFT = IOFFT + DISSIZ * NSUM
  250 CONTINUE
C
C     Summation into expanded D3T3 : D3T3EXP(A,Bc) = W(A,E) * T3EXP(Bc,E)
C                                    D3T3EXP(a,bC) = W(a,e) * T3EXP(bC,e)
C
      LISTW = 57 + ISPIN1
      CALL GETLST(CORE(I000),MK,1,2,IRPMK,LISTW)
C
      CALL SYMTRW2(CORE(I010),CORE(I020),CORE(I030),IADT,IADTX,IRPIJM,
     1             ISPIN1,ISPIN2)
C
      IOFFT = I020
      DO 260 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN1).EQ.0) GOTO 260
C
      IRPA  = DIRPRD(IRPE,IRPMK)
      IRPAE = IRPMK
      IOFFW = IOFFVV(IRPA,IRPAE,ISPIN1 + 2) + 1
      IRPBC = DIRPRD(IRPE,IRPIJM)
      DISSIZ = VRT(IRPA,ISPIN1)
      NDIS   = IRPDPD(IRPBC,13)
      NSUM   = VRT(IRPE,ISPIN1)
      IOFFDT = IADDTX(IRPBC)
C
      CALL XGEMM('T','T',DISSIZ,NDIS,NSUM,SIGN,
     1           CORE(IOFFW),NSUM,CORE(IOFFT),NDIS,1.0D+00,
     1           D3T3EXP(IOFFDT),DISSIZ)
C
      IOFFT = IOFFT + NDIS * NSUM
  260 CONTINUE
C
  290 CONTINUE
  300 CONTINUE
C
C     Contributions to case 2 from case 3 and vice versa.
C
C     I000 - contribution to D3T3EXP.
C
      LENINC = 0
      DO 305 IRPBC=1,NIRREP
      IRPA = DIRPRD(IRPBC,IRPIJK)
      LENINC = LENINC + VRT(IRPA,ISPIN1) * IRPDPD(IRPBC,13)
  305 CONTINUE
      CALL ZERO(CORE(1),LENINC)
C
      DO 400 IRPM=1,NIRREP
C
      IF(POP(IRPM,ISPIN2).EQ.0)                  GOTO 400
      IF(POP(IRPM,ISPIN2).LT.2.AND.IRPM.EQ.IRPK) GOTO 400
C
      IF(IRPM.LE.IRPK) MKJPOS = IJKPOS(IRPM,IRPK,IRPJ,2)
      IF(IRPM.GT.IRPK) MKJPOS = IJKPOS(IRPK,IRPM,IRPJ,2)
C
      IRPKJ  = DIRPRD(IRPK,IRPJ)
      IRPMKJ = DIRPRD(IRPM,IRPKJ)
C
      IOFF = IJKOFF(MKJPOS,IRPMKJ,ISPIN2 + 1)
C
      LN    = DISTSZ(IRPMKJ,ISPIN2 + 1 + 4)
      LNEXP = 0
      DO 310 IRREP=1,NIRREP
      JRREP = DIRPRD(IRREP,IRPMKJ)
      LNEXP = LNEXP + IRPDPD(JRREP,13) * VRT(IRREP,ISPIN2)
      LENT(IRREP) = IRPDPD(JRREP,ISPIN2) * VRT(IRREP,ISPIN1)
  310 CONTINUE
C
      DO 320 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADT(IRREP) = 1
      ELSE
      IADT(IRREP) = IADT(IRREP-1) + LENT(IRREP-1)
      ENDIF
  320 CONTINUE
C
C     I000 -  W(C,E)
C     I010 - T3(A<B,e)
C     I020 - T3(B,e,A) (labelled by A)
C
      IRPMI = DIRPRD(IRPI,IRPM)
      I000 = 1
      I010 = I000 + LENINC
      I020 = I010 + IRPDPD(IRPMI,13)
      I030 = I020 + LN
      I040 = I030 + LNEXP
      I050 = I040 + LNEXP
      NEED = IINTFP * I050
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
      CALL INSMEM('T3WT323R',NEED,MAXCOR)
      ENDIF
C
      DO  390    M=1,POP(IRPM,ISPIN2)
C
      IF(IRPM.EQ.IRPK.AND.M.EQ.K) GOTO 390
C
      SIGN =  1.0D+00
C
      IF(IRPM.EQ.IRPK)THEN
      NMK = (POP(IRPM,ISPIN2) * (POP(IRPM,ISPIN2) - 1))/2
       IF(M.LT.K)THEN
       MK = INDEX(K-1) + M
       ELSE
       MK = INDEX(M-1) + K
       SIGN = - SIGN
       ENDIF
      ELSE
      NMK = POP(IRPM,ISPIN2) * POP(IRPK,ISPIN2)
       IF(IRPM.LT.IRPK)THEN
       MK = (K-1)*POP(IRPM,ISPIN2) + M
       ELSE
       MK = (M-1)*POP(IRPK,ISPIN2) + K
       SIGN = - SIGN
       ENDIF
      ENDIF
      MKJVAL = IOFF + (J-1)*NMK + MK
C
      IF(IUHF.EQ.0)THEN
      CALL GETLIST(CORE(I020),MKJVAL,1,1,IRPMKJ,ISPIN1 + 1 + 4)
      ELSE
      CALL GETLIST(CORE(I020),MKJVAL,1,1,IRPMKJ,ISPIN2 + 1 + 4)
      ENDIF
      CALL SYMTRW2(CORE(I020),CORE(I030),CORE(I040),IADT,IADTX,IRPMKJ,
     1             ISPIN2,ISPIN1)
C      write(6,*) ' t3, t3exp, t3exp ',M,I,J,K,MKJVAL
C      CALL SUMBLK(CORE(I020),LN)
C      CALL SUMBLK(CORE(I030),LNEXP)
C      CALL SUMBLK(CORE(I040),LNEXP)
C
C     We have (cB,e) or (Cb,E) at I020. We need (Bc,e) or (bC,E).
C     The transposition is done for each irpe below.
C
      MI = IOFFOO(IRPI,IRPMI,4+ISPIN2) + (I-1)*POP(IRPM,ISPIN2) + M
      IF(IUHF.EQ.0)THEN
      LISTW = 56
      ELSE
      LISTW = 58 - ISPIN1
      ENDIF
      CALL GETLST(CORE(I010),MI,1,2,IRPMI,LISTW)
C      write(6,*) ' intermediate '
C      CALL SUMBLK(CORE(I010),IRPDPD(IRPMI,13))
C
C     Summation into expanded D3T3 : D3T3EXP(A,Bc) = W(A,E) * T3EXP(Bc,E)
C                                    D3T3EXP(a,bC) = W(a,e) * T3EXP(bC,e)
C
C
      IOFFT = I030
      DO 360 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN2).EQ.0) GOTO 360
C
      IRPA  = DIRPRD(IRPE,IRPMI)
      IRPAE = IRPMI
      IOFFW = IOFFVV(IRPA,IRPAE,ISPIN2 + 4) + I010
      IRPBC = DIRPRD(IRPE,IRPMKJ)
      DISSIZ = VRT(IRPA,ISPIN1)
      NDIS   = IRPDPD(IRPBC,13)
      NSUM   = VRT(IRPE,ISPIN2)
      IOFFDT = IADDTX(IRPBC)
C
C     Do not do this transposition now. instead, form (A,cB) and transpose
C     this after all the summations have been performed.
c      I040 = I030 + IRPDPD(IRPBC,13)
c      I050 = I040 + IRPDPD(IRPBC,13)
c      I060 = I050 + IRPDPD(IRPBC,13)
c      CALL SYMTR3(IRPBC,VRT(1,ISPIN2),VRT(1,ISPIN1),IRPDPD(IRPBC,13),
c     1            VRT(IRPE,ISPIN2),CORE(IOFFT),CORE(I030),CORE(I040),
c     1            CORE(I050))
C
c      CALL XGEMM('T','T',DISSIZ,NDIS,NSUM,SIGN,
c     1           CORE(IOFFW),NSUM,CORE(IOFFT),NDIS,1.0D+00,
c     1           D3T3EXP(IOFFDT),DISSIZ)
      CALL XGEMM('T','T',DISSIZ,NDIS,NSUM,SIGN,
     1           CORE(IOFFW),NSUM,CORE(IOFFT),NDIS,1.0D+00,
     1           CORE(IOFFDT),DISSIZ)
C
      IOFFT = IOFFT + NDIS * NSUM
  360 CONTINUE
  390 CONTINUE
  400 CONTINUE
C
c      DO  405 IRPA=1,NIRREP
c      IRPBC  = DIRPRD(IRPA,IRPIJK)
c      DISSIZ = VRT(IRPA,ISPIN1)
c      NDIS   = IRPDPD(IRPBC,13)
cC      
c      I020 = I010 + MAX(DISSIZ,NDIS)
c      I030 = I020 + MAX(DISSIZ,NDIS)
c      CALL SYMTR1(IRPBC,VRT(1,ISPIN2),VRT(1,ISPIN1),DISSIZ,
c     1            CORE(IADDTX(IRPBC)),
c     1            CORE(I010),CORE(I020),CORE(I030))
cC
c      CALL VADD(D3T3EXP(IADDTX(IRPBC)),D3T3EXP(IADDTX(IRPBC)),
c     1          CORE(IADDTX(IRPBC)),DISSIZ*NDIS, 1.0D+00)
c  405 CONTINUE
C
      DO 500 IRPM=1,NIRREP
C
      IF(POP(IRPM,ISPIN2).EQ.0)                  GOTO 500
      IF(POP(IRPM,ISPIN2).LT.2.AND.IRPM.EQ.IRPK) GOTO 500
C
      IF(IRPM.LE.IRPK) MKIPOS = IJKPOS(IRPM,IRPK,IRPI,2)
      IF(IRPM.GT.IRPK) MKIPOS = IJKPOS(IRPK,IRPM,IRPI,2)
C
      IRPKI  = DIRPRD(IRPK,IRPI)
      IRPMKI = DIRPRD(IRPM,IRPKI)
C
      IOFF = IJKOFF(MKIPOS,IRPMKI,ISPIN2 + 1)
C
      LN    = DISTSZ(IRPMKI,ISPIN2 + 1 + 4)
      LNEXP = 0
      DO 410 IRREP=1,NIRREP
      JRREP = DIRPRD(IRREP,IRPMKI)
      LNEXP = LNEXP + IRPDPD(JRREP,13) * VRT(IRREP,ISPIN2)
      LENT(IRREP) = IRPDPD(JRREP,ISPIN2) * VRT(IRREP,ISPIN1)
  410 CONTINUE
C
      DO 420 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADT(IRREP) = 1
      ELSE
      IADT(IRREP) = IADT(IRREP-1) + LENT(IRREP-1)
      ENDIF
  420 CONTINUE
C
C     I000 -  W(C,E)
C     I010 - T3(A<B,e)
C     I020 - T3(B,e,A) (labelled by A)
C
      IRPMJ = DIRPRD(IRPJ,IRPM)
      I000 = 1
      I010 = I000 + LENINC
      I020 = I010 + IRPDPD(IRPMJ,13)
      I030 = I020 + LN
      I040 = I030 + LNEXP
      I050 = I040 + LNEXP
      NEED = IINTFP * I050
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
      CALL INSMEM('T3WT323R',NEED,MAXCOR)
      ENDIF
C
      DO  490    M=1,POP(IRPM,ISPIN2)
C
      IF(IRPM.EQ.IRPK.AND.M.EQ.K) GOTO 490
C
      SIGN = -1.0D+00
C
      IF(IRPM.EQ.IRPK)THEN
      NMK = (POP(IRPM,ISPIN2) * (POP(IRPM,ISPIN2) - 1))/2
       IF(M.LT.K)THEN
       MK = INDEX(K-1) + M
       ELSE
       MK = INDEX(M-1) + K
       SIGN = - SIGN
       ENDIF
      ELSE
      NMK = POP(IRPM,ISPIN2) * POP(IRPK,ISPIN2)
       IF(IRPM.LT.IRPK)THEN
       MK = (K-1)*POP(IRPM,ISPIN2) + M
       ELSE
       MK = (M-1)*POP(IRPK,ISPIN2) + K
       SIGN = - SIGN
       ENDIF
      ENDIF
      MKIVAL = IOFF + (I-1)*NMK + MK
C
      IF(IUHF.EQ.0)THEN
      CALL GETLIST(CORE(I020),MKIVAL,1,1,IRPMKI,ISPIN1 + 1 + 4)
      ELSE
      CALL GETLIST(CORE(I020),MKIVAL,1,1,IRPMKI,ISPIN2 + 1 + 4)
      ENDIF
      CALL SYMTRW2(CORE(I020),CORE(I030),CORE(I040),IADT,IADTX,IRPMKI,
     1             ISPIN2,ISPIN1)
C
C     We have (cB,e) or (Cb,E) at I020. We need (Bc,e) or (bC,E).
C     The transposition is done for each irpe below.
C
      MJ = IOFFOO(IRPJ,IRPMJ,4+ISPIN2) + (J-1)*POP(IRPM,ISPIN2) + M
      IF(IUHF.EQ.0)THEN
      LISTW = 56
      ELSE
      LISTW = 58 - ISPIN1
      ENDIF
      CALL GETLST(CORE(I010),MJ,1,2,IRPMJ,LISTW)
C
C     Summation into expanded D3T3 : D3T3EXP(A,Bc) = W(A,E) * T3EXP(Bc,E)
C                                    D3T3EXP(a,bC) = W(a,e) * T3EXP(bC,e)
C
C
      IOFFT = I030
      DO 460 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN2).EQ.0) GOTO 460
C
      IRPA  = DIRPRD(IRPE,IRPMJ)
      IRPAE = IRPMJ
      IOFFW = IOFFVV(IRPA,IRPAE,ISPIN2 + 4) + I010
      IRPBC = DIRPRD(IRPE,IRPMKI)
      DISSIZ = VRT(IRPA,ISPIN1)
      NDIS   = IRPDPD(IRPBC,13)
      NSUM   = VRT(IRPE,ISPIN2)
      IOFFDT = IADDTX(IRPBC)
C
c      I040 = I030 + IRPDPD(IRPBC,13)
c      I050 = I040 + IRPDPD(IRPBC,13)
c      I060 = I050 + IRPDPD(IRPBC,13)
c      CALL SYMTR3(IRPBC,VRT(1,ISPIN2),VRT(1,ISPIN1),IRPDPD(IRPBC,13),
c     1            VRT(IRPE,ISPIN2),CORE(IOFFT),CORE(I030),CORE(I040),
c     1            CORE(I050))
C
c      CALL XGEMM('T','T',DISSIZ,NDIS,NSUM,SIGN,
c     1           CORE(IOFFW),NSUM,CORE(IOFFT),NDIS,1.0D+00,
c     1           D3T3EXP(IOFFDT),DISSIZ)
      CALL XGEMM('T','T',DISSIZ,NDIS,NSUM,SIGN,
     1           CORE(IOFFW),NSUM,CORE(IOFFT),NDIS,1.0D+00,
     1           CORE(IOFFDT),DISSIZ)
C
      IOFFT = IOFFT + NDIS * NSUM
  460 CONTINUE
  490 CONTINUE
  500 CONTINUE
C
      DO  405 IRPA=1,NIRREP
      IRPBC  = DIRPRD(IRPA,IRPIJK)
      DISSIZ = VRT(IRPA,ISPIN1)
      NDIS   = IRPDPD(IRPBC,13)
C      
      I020 = I010 + MAX(DISSIZ,NDIS)
      I030 = I020 + MAX(DISSIZ,NDIS)
      CALL SYMTR1(IRPBC,VRT(1,ISPIN2),VRT(1,ISPIN1),DISSIZ,
     1            CORE(IADDTX(IRPBC)),
     1            CORE(I010),CORE(I020),CORE(I030))
C
      CALL VADD(D3T3EXP(IADDTX(IRPBC)),D3T3EXP(IADDTX(IRPBC)),
     1          CORE(IADDTX(IRPBC)),DISSIZ*NDIS, 1.0D+00)
  405 CONTINUE
C
C     Case 1/4 contribution to 2/3.
C
C     ABc IJk  =  ABE IJM   *  Mc Ek   (D3T3(A<B,c) = T3EXP(A<B,E)*W(E,c))
C     abC ijK  =  abe ijm   *  mC eK   (D3T3(a<b,C) = T3EXP(a<b,e)*W(e,C))
C
      DO 600 IRPM=1,NIRREP
C
      IF(POP(IRPM,ISPIN1).EQ.0) GOTO 600
      IF(POP(IRPM,ISPIN1).LT.2.AND.(IRPM.EQ.IRPI.OR .IRPM.EQ.IRPJ))
     1                          GOTO 600
      IF(POP(IRPM,ISPIN1).LT.3.AND. IRPM.EQ.IRPI.AND.IRPM.EQ.IRPJ )
     1                          GOTO 600
C
      IF(IRPM.GE.IRPJ) IJMPOS = IJKPOS(IRPI,IRPJ,IRPM,1)
      IF(IRPM.LE.IRPI) IJMPOS = IJKPOS(IRPM,IRPI,IRPJ,1)
      IF(IRPM.GT.IRPI.AND.IRPM.LT.IRPJ)
     1                 IJMPOS = IJKPOS(IRPI,IRPM,IRPJ,1)
C
      IRPIJ  = DIRPRD(IRPI,IRPJ)
      IRPIJM = DIRPRD(IRPM,IRPIJ)
C
      IOFF = IJKOFF(IJMPOS,IRPIJM,ISPIN1 + 2*(ISPIN1-1))
C
      LN    = DISTSZ(IRPIJM,ISPIN1 + 2*(ISPIN1-1))
      N14   = NDISTS(IRPIJM,ISPIN1 + 2*(ISPIN1-1))
C
C     Skip if there are no AAA/BBB amplitudes.
C
      IF(LN.EQ.0.OR.N14.EQ.0) GOTO 600
C
      LNEXP = 0
      DO 510 IRREP=1,NIRREP
      JRREP = DIRPRD(IRREP,IRPIJM)
      LNEXP = LNEXP + IRPDPD(JRREP,ISPIN1) * VRT(IRREP,ISPIN1)
      LENT(IRREP)   = IRPDPD(JRREP,ISPIN1) * VRT(IRREP,ISPIN1)
  510 CONTINUE
C
      DO 520 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADT(IRREP) = 1
      ELSE
      IADT(IRREP) = IADT(IRREP-1) + LENT(IRREP-1)
      ENDIF
  520 CONTINUE
C
C     I000 -  W(C,E)
C     I010 - T3(A<B<E)
C     I020 - T3(A<B,E) (labelled by E)
C
      IRPMK = DIRPRD(IRPK,IRPM)
      I000 = 1
      I010 = I000 + IRPDPD(IRPMK,13)
      I020 = I010 + LN
      I030 = I020 + LNEXP
      I040 = I030 + LNEXP
      NEED = IINTFP * I040
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
      CALL INSMEM('T3WT323R',NEED,MAXCOR)
      ENDIF
C
      NONEQL = .FALSE.
      IJMEQL = .FALSE.
      IJEQL  = .FALSE.
      IMEQL  = .FALSE.
      JMEQL  = .FALSE.
      NONEQL = IRPI.NE.IRPJ.AND.IRPI.NE.IRPM.AND.IRPJ.NE.IRPM
      IJMEQL = IRPI.EQ.IRPJ.AND.IRPI.EQ.IRPM
      IJEQL  = IRPI.EQ.IRPJ.AND.IRPI.NE.IRPM
      IMEQL  = IRPI.EQ.IRPM.AND.IRPI.NE.IRPJ
      JMEQL  = IRPJ.EQ.IRPM.AND.IRPJ.NE.IRPI
C
      DO  590    M=1,POP(IRPM,ISPIN1)
C
      IF(IRPM.EQ.IRPI.AND.M.EQ.I) GOTO 590
      IF(IRPM.EQ.IRPJ.AND.M.EQ.J) GOTO 590
C
      SIGN = 1.0D+00
C
      IF(NONEQL)THEN
       IF(IRPM.LT.IRPI)THEN
       IJMVAL = IOFF + (J-1)*POP(IRPI,ISPIN1)*POP(IRPM,ISPIN1) +
     1                 (I-1)*POP(IRPM,ISPIN1) + M
       ENDIF
       IF(IRPM.GT.IRPJ)THEN
       IJMVAL = IOFF + (M-1)*POP(IRPJ,ISPIN1)*POP(IRPI,ISPIN1) +
     1                 (J-1)*POP(IRPI,ISPIN1) + I
       ENDIF
       IF(IRPM.GT.IRPI.AND.IRPM.LT.IRPJ)THEN
       IJMVAL = IOFF + (J-1)*POP(IRPM,ISPIN1)*POP(IRPI,ISPIN1) +
     1                 (M-1)*POP(IRPI,ISPIN1) + I
       SIGN = - SIGN
       ENDIF
      ENDIF
C
      IF(IJMEQL)THEN
      MAXVAL = MAX(M,I,J)
      MINVAL = MIN(M,I,J)
      IF(M.GT.J) MIDVAL = J
      IF(M.LT.I) MIDVAL = I
      IF(M.GT.I.AND.M.LT.J)THEN
      MIDVAL = M
      SIGN = - SIGN
      ENDIF
      IJMVAL = IOFF + ((MAXVAL-1)*(MAXVAL-2)*(MAXVAL-3))/6 +
     1                INDEX(MIDVAL-1) + MINVAL
      ENDIF
C
      IF(IMEQL)THEN
       IF(M.LT.I)THEN
       MI = INDEX(I-1) + M
       ELSE
       MI = INDEX(M-1) + I
       SIGN = - SIGN
       ENDIF
       IJMVAL = IOFF + (J-1)*((POP(IRPM,ISPIN1)*(POP(IRPM,ISPIN1)-1))/2)
     1               + MI
      ENDIF
C
      IF(JMEQL)THEN
       IF(M.LT.J)THEN
       MJ = INDEX(J-1) + M
       SIGN = - SIGN
       ELSE
       MJ = INDEX(M-1) + J
       ENDIF
       IJMVAL = IOFF + (MJ-1)*POP(IRPI,ISPIN1) + I
      ENDIF
C
      IF(IJEQL)THEN
       IF(IRPM.GT.IRPJ)THEN
       IJMVAL = IOFF + (M-1)*((POP(IRPJ,ISPIN1)*(POP(IRPJ,ISPIN1)-1))/2)
     1               + INDEX(J-1) + I
       ELSE
       IJMVAL = IOFF + (INDEX(J-1) + I - 1)*POP(IRPM,ISPIN1) + M
       ENDIF
      ENDIF
C
      CALL GETLIST(CORE(I010),IJMVAL,1,1,IRPIJM,1 + 3*(ISPIN1-1) + 4)
C
      CALL ZERO(CORE(I020),LNEXP)
      CALL SYMEXPT3(CORE(I010),CORE(I020),IADT,ISPIN1,IRPIJM)
C
      MK = IOFFOO(IRPK,IRPMK,4+ISPIN1) + (K-1)*POP(IRPM,ISPIN1) + M
      LISTW = 55 + ISPIN1
      CALL GETLST(CORE(I000),MK,1,2,IRPMK,LISTW)
C
      IOFFT = I020
      DO 560 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN1).EQ.0) GOTO 560
C
      IRPC  = DIRPRD(IRPE,IRPMK)
      IRPCE = IRPMK
      IOFFW = IOFFVV(IRPC,IRPCE,ISPIN1 + 4) + 1
      IRPAB = DIRPRD(IRPE,IRPIJM)
      DISSIZ = IRPDPD(IRPAB,ISPIN1)
      NDIS   = VRT(IRPC,ISPIN2)
      NSUM   = VRT(IRPE,ISPIN1)
      IOFFDT = IADDT(IRPC)
C
      CALL XGEMM('N','N',DISSIZ,NDIS,NSUM,SIGN,
     1           CORE(IOFFT),DISSIZ,CORE(IOFFW),NSUM,1.0D+00,
     1           D3T3(IOFFDT),DISSIZ)
C
      IOFFT = IOFFT + DISSIZ * NSUM
  560 CONTINUE
  590 CONTINUE
  600 CONTINUE
      RETURN
      END