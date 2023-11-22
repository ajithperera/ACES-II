      SUBROUTINE GT3WT23(T3,W,CORE,
     &                   IADT3,IADW,
     &                   I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPJK,IRPIK,
     &                   IRREPT,LIST2OFF,IRREPW,LWOOFF,LWVOFF,
     &                   MAXCOR)
C
C-----------------------------------------------------------------------
C
C     This is moving toward being a generic subroutine for computing
C     contractions of the form Z3 = W*T2 for the BBA spin case of Z3.
C
C     T3    Z3(a<b,C). Labelled by IRPC. Must be initialized by caller.
C     W     Z3(a,bC).  Labelled by IRPBC. Must be initialized by caller.
C           [N.B. W is not the "W" which is contracted with T2]
C     CORE  Scratch space.
C     IADT3 Addresses for symmetry blocks of T3. Must be supplied by
C           caller.
C     IADW  Addresses for symmetry blocks of W. Must be supplied by
C           caller.
C     
C     List expectations (Watch for diversions from standards)
C
C     1-NIRREP  LIST2OFF + 1     T2(AB,IJ)      A<B ; I<J
C     1-NIRREP  LIST2OFF + 2     T2(ab,ij)      a<b ; i<j
C     1-NIRREP  LIST2OFF + 3     T2(Ab,Ij)      b,A ; I;j
C
C     1-NIRREP  LWOOFF + 1        W(KA,IJ)      K,A ; I<J
C     1-NIRREP  LWOOFF + 2        W(ka,ij)      k,a ; i<j
C     1-NIRREP  LWOOFF + 3        W(Ak,Ij)      A,k ; I,j
C     1-NIRREP  LWOOFF + 4        W(Ka,Ij)      K,a ; I,j
C
C     1-NIRREP  LWVOFF + 1        W(AB,CI)      A<B ; C,I
C     1-NIRREP  LWVOFF + 2        W(ab,ci)      a<b ; c,i
C     1-NIRREP  LWVOFF + 3        W(Ab,Ic)      b,A ; c,I
C     1-NIRREP  LWVOFF + 4        W(Ab,Ci)      b,A ; C,i
C
C     The RHF/UHF question is moot for this subroutine. It should only
C     be called in UHF calculations, as we only generate AAB triples in
C     RHF. Hence this subroutine assumes the full UHF set of lists is
C     available.
C
C     W (Hbar elements, integrals) are assumed to be totally symmetric.
C     T3 and T2 have overall symmetry IRREPX.
C-----------------------------------------------------------------------
C
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE
C-----------------------------------------------------------------------
      DOUBLE PRECISION T3,W,CORE
      INTEGER IADT3,IADW,I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPJK,IRPIK,
     &        IRREPT,LIST2OFF,IRREPW,LWOOFF,LWVOFF,MAXCOR
C-----------------------------------------------------------------------
      DOUBLE PRECISION ONE,ONEM
      INTEGER IINTLN,IFLTLN,IINTFP,IALONE,IBITWD,
     &        NSTART,NIRREP,IRREPA,IRREPB,DIRPRD,IRPDPD,ISYTYP,ID,
     &        POP,VRT,NTAA,NTBB,NF1AA,NF2AA,NF1BB,NF2BB,
     &        IOFFOO,IOFFVV,IOFFVO,INDEX
C-----------------------------------------------------------------------
      INTEGER IADT2,LENT2,DSZ,DSZEXP,
     &        IRPA,IRPAB,IRPAE,IRPBC,IRPC,IRPE,IRPEC,IRPEI,IRPEJ,IRPEK,
     &        IRPM,IRPMA,IRPMC,IRPMI,IRPMJ,IRPMK,IRPIM,IRPJM,
     &        IJ,KJ,KI,M,IMOFF,JMOFF,JM,IM,
     &        I000,I010,I020,I030,I040,I050,NEED,LENV,IOFFV,
     &        ISPIN1,ISPIN2,
     &        LISWOA,LISWOB,LISWOC,LISWOD,LISWVA,LISWVB,LISWVC,LISWVD
C-----------------------------------------------------------------------
      DIMENSION T3(1),W(1),CORE(1)
      DIMENSION IADT3(8),IADW(8) 
C
      DIMENSION IADT2(8) ,LENT2(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     &                NF1BB,NF2BB
C
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
      DATA ONE /1.0D+00/, ONEM /-1.0D+00/
C
      INDEX(I) = I*(I-1)/2
C
      LISWVA = LWVOFF + 1
      LISWVB = LWVOFF + 2
      LISWVC = LWVOFF + 3
      LISWVD = LWVOFF + 4
C
      LISWOA = LWOOFF + 1
      LISWOB = LWOOFF + 2
      LISWOC = LWOOFF + 3
      LISWOD = LWOOFF + 4
C
      ISPIN1 = 2
      ISPIN2 = 1
C
C-----------------------------------------------------------------------
C      ae                be
C     T   <Cb//Ke>  -   T   <Ca//Ke>
C      ij                ij
C
C-----------------------------------------------------------------------
C
      IRPAE = DIRPRD(IRREPT,IRPIJ)
C     
      DSZ    = IRPDPD(IRPAE,     ISPIN1)
      DSZEXP = IRPDPD(IRPAE,18 + ISPIN1)
C   
      I000 = 1
      I010 = I000 + DSZEXP
      I020 = I010 + DSZ
C
      IF(IRPI.EQ.IRPJ)THEN
      IJ = IOFFOO(IRPJ,IRPIJ,ISPIN1) + INDEX(J-1) + I
      ELSE
      IJ = IOFFOO(IRPJ,IRPIJ,ISPIN1) + (J-1)*POP(IRPI,ISPIN1) + I
      ENDIF
      CALL GETLST(CORE(I010),IJ,1,1,IRPIJ,LIST2OFF + ISPIN1)
      CALL SYMEXP2(IRPAE,VRT(1,ISPIN1),DSZEXP,DSZ,1,CORE(I000),
     &                                              CORE(I010))
C
      DO   10 IRPE=1,NIRREP
        IRPA = DIRPRD(IRPE,IRPAE)
        LENT2(IRPE) = VRT(IRPE,ISPIN1) * VRT(IRPA,ISPIN1)
          IF(IRPE.EQ.1)THEN
            IADT2(IRPE) = 1
          ELSE
            IADT2(IRPE) = IADT2(IRPE-1) + LENT2(IRPE-1)
          ENDIF
   10 CONTINUE
C
      DO   20 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN1).EQ.0) GOTO 20

      IRPA  = DIRPRD(IRPE,IRPAE)
      IRPEK = DIRPRD(IRPE,IRPK)
      IRPBC = DIRPRD(IRREPW,IRPEK)
      LENV  = IRPDPD(IRPBC,13) * VRT(IRPE,ISPIN1)
      I020  = I010 + LENV
      NEED  = IINTFP * I020
C
      IF(NEED.GT.MAXCOR)THEN
        WRITE(6,1010)
        CALL INSMEM('GT3WT23',NEED,MAXCOR)
      ENDIF
C
      CALL GETLST(CORE(I010),
     &            IOFFVO(IRPK,IRPEK,4)+(K-1)*VRT(IRPE,ISPIN1)+1,
     &            VRT(IRPE,ISPIN1),2,IRPEK,LISWVC)
C
      CALL XGEMM('N','T',
     &           VRT(IRPA,ISPIN1),IRPDPD(IRPBC,13),VRT(IRPE,ISPIN1),
     &           ONE,CORE(I000 + IADT2(IRPE) - 1),VRT(IRPA,ISPIN1),
     &               CORE(I010),IRPDPD(IRPBC,13),1.0D+00,
     &           W(IADW(IRPBC)),VRT(IRPA,ISPIN1))
C
   20 CONTINUE
C
C-----------------------------------------------------------------------
C
C     - T(Ea,Kj) * W(Cb,Ei)   +  T(Eb,Kj) * W(Ca,Ei)
C
C-----------------------------------------------------------------------
C
      IRPAE = DIRPRD(IRREPT,IRPJK)
C
      DO  110 IRPE=1,NIRREP
        IRPA = DIRPRD(IRPE,IRPAE)
        LENT2(IRPE) = VRT(IRPE,ISPIN2) * VRT(IRPA,ISPIN1)
          IF(IRPE.EQ.1)THEN
            IADT2(IRPE) = 1
          ELSE
            IADT2(IRPE) = IADT2(IRPE-1) + LENT2(IRPE-1)
          ENDIF
  110 CONTINUE
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPAE,13)
      KJ = IOFFOO(IRPJ,IRPJK,5) + (J-1)*POP(IRPK,ISPIN2) + K
      CALL GETLST(CORE(I000),KJ,1,1,IRPJK,LIST2OFF + 3)
C
      DO  120 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN2).EQ.0) GOTO 120
C
      IRPA  = DIRPRD(IRPE,IRPAE)
      IRPEI = DIRPRD(IRPE,IRPI)
      IRPBC = DIRPRD(IRREPW,IRPEI)
      LENV  = IRPDPD(IRPBC,13) * VRT(IRPE,ISPIN2)
C
      I020 = I010 + IRPDPD(IRPBC,13)
      I030 = I020 + IRPDPD(IRPBC,13)
      I040 = I030 + IRPDPD(IRPBC,13)
      I050 = I040 + LENV
      NEED = IINTFP * I050
C
      IF(NEED.GT.MAXCOR)THEN
        WRITE(6,1010)
        CALL INSMEM('GT3WT23',NEED,MAXCOR)
      ENDIF
C
      CALL GETLST(CORE(I040),
     &            IOFFVO(IRPI,IRPEI,4) + (I-1)*VRT(IRPE,ISPIN2) + 1,
     &            VRT(IRPE,ISPIN2),2,IRPEI,LISWVD)
C
      CALL XGEMM('N','T',
     &           VRT(IRPA,ISPIN1),IRPDPD(IRPBC,13),VRT(IRPE,ISPIN2),
     &           ONEM,
     &           CORE(I000 + IADT2(IRPE) - 1),VRT(IRPA,ISPIN1),
     &           CORE(I040),IRPDPD(IRPBC,13),ONE,
     &           W(IADW(IRPBC)),VRT(IRPA,ISPIN1))
  120 CONTINUE
C
C-----------------------------------------------------------------------
C
C       T(Ea,Ki) * W(Cb,Ej)   -  T(Eb,Ki) * W(Ca,Ej)
C
C-----------------------------------------------------------------------
C
      IRPAE = DIRPRD(IRREPT,IRPIK)
C
      DO  210 IRPE=1,NIRREP
        IRPA = DIRPRD(IRPE,IRPAE)
        LENT2(IRPE) = VRT(IRPE,ISPIN2) * VRT(IRPA,ISPIN1)
          IF(IRPE.EQ.1)THEN
            IADT2(IRPE) = 1
          ELSE
            IADT2(IRPE) = IADT2(IRPE-1) + LENT2(IRPE-1)
          ENDIF
  210 CONTINUE
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPAE,13)
      KI = IOFFOO(IRPI,IRPIK,5) + (I-1)*POP(IRPK,ISPIN2) + K
      CALL GETLST(CORE(I000),KI,1,1,IRPIK,LIST2OFF + 3)
C
      DO  220 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN2).EQ.0) GOTO 220
C
      IRPA  = DIRPRD(IRPE,IRPAE)
      IRPEJ = DIRPRD(IRPE,IRPJ)
      IRPBC = DIRPRD(IRREPW,IRPEJ)
      LENV  = IRPDPD(IRPBC,13) * VRT(IRPE,ISPIN2)
C
      I020 = I010 + IRPDPD(IRPBC,13)
      I030 = I020 + IRPDPD(IRPBC,13)
      I040 = I030 + IRPDPD(IRPBC,13)
      I050 = I040 + LENV
      NEED = IINTFP * I050
C
      IF(NEED.GT.MAXCOR)THEN
        WRITE(6,1010)
        CALL INSMEM('GT3WT23',NEED,MAXCOR)
      ENDIF
C
      CALL GETLST(CORE(I040),
     &            IOFFVO(IRPJ,IRPEJ,4) + (J-1)*VRT(IRPE,ISPIN2) + 1,
     &            VRT(IRPE,ISPIN2),2,IRPEJ,LISWVD)
C
      CALL XGEMM('N','T',
     &           VRT(IRPA,ISPIN1),IRPDPD(IRPBC,13),VRT(IRPE,ISPIN2),
     &           ONE,
     &           CORE(I000 + IADT2(IRPE) - 1),VRT(IRPA,ISPIN1),
     &           CORE(I040),IRPDPD(IRPBC,13),ONE,
     &           W(IADW(IRPBC)),VRT(IRPA,ISPIN1))
  220 CONTINUE
C
C-----------------------------------------------------------------------
C
C     T(Ce,Ki) * W(ab,ej)
C
C-----------------------------------------------------------------------
C
      IRPEC = DIRPRD(IRREPT,IRPIK)
C
      DO  310 IRPC=1,NIRREP
        IRPE = DIRPRD(IRPEC,IRPC)
        LENT2(IRPC) = VRT(IRPC,ISPIN2) * VRT(IRPE,ISPIN1)
          IF(IRPC.EQ.1)THEN
            IADT2(IRPC) = 1
          ELSE
            IADT2(IRPC) = IADT2(IRPC-1) + LENT2(IRPC-1)
          ENDIF
  310 CONTINUE
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPEC,13)
      KI = IOFFOO(IRPI,IRPIK,5) + (I-1)*POP(IRPK,ISPIN2) + K
      CALL GETLST(CORE(I000),KI,1,1,IRPIK,LIST2OFF + 3)
C
      DO  320 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN1).EQ.0) GOTO 320
C
      IRPC  = DIRPRD(IRPE,IRPEC)
      IRPEJ = DIRPRD(IRPE,IRPJ)
      IRPAB = DIRPRD(IRREPW,IRPEJ)
C
      I020 = I010 + IRPDPD(IRPAB,ISPIN1) * VRT(IRPE,ISPIN1)
      NEED = IINTFP * I020
C
      IF(NEED.GT.MAXCOR)THEN
        WRITE(6,1010)
        CALL INSMEM('GT3WT23',NEED,MAXCOR)
      ENDIF
C
      CALL GETLST(CORE(I010),
     &       IOFFVO(IRPJ,IRPEJ,ISPIN1) + (J-1)*VRT(IRPE,ISPIN1) + 1,
     &            VRT(IRPE,ISPIN1),2,IRPEJ,LISWVB)
C
      CALL XGEMM('N','N',
     &           IRPDPD(IRPAB,ISPIN1),VRT(IRPC,ISPIN2),VRT(IRPE,ISPIN1),
     &           ONE,
     &           CORE(I010),IRPDPD(IRPAB,ISPIN1),
     &           CORE(I000 + IADT2(IRPC) - 1),VRT(IRPE,ISPIN1),
     &           ONE,
     &           T3(IADT3(IRPC)),IRPDPD(IRPAB,ISPIN1))
C  
  320 CONTINUE
C
C-----------------------------------------------------------------------
C
C   - T(Ce,Kj) * W(ab,ei)
C
C-----------------------------------------------------------------------
C
      IRPEC = DIRPRD(IRREPT,IRPJK)
C
      DO  410 IRPC=1,NIRREP
        IRPE = DIRPRD(IRPEC,IRPC)
        LENT2(IRPC) = VRT(IRPC,ISPIN2) * VRT(IRPE,ISPIN1)
          IF(IRPC.EQ.1)THEN
            IADT2(IRPC) = 1
          ELSE
            IADT2(IRPC) = IADT2(IRPC-1) + LENT2(IRPC-1)
          ENDIF
  410 CONTINUE
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPEC,13)
      KJ = IOFFOO(IRPJ,IRPJK,5) + (J-1)*POP(IRPK,ISPIN2) + K
      CALL GETLST(CORE(I000),KJ,1,1,IRPJK,LIST2OFF + 3)
C
      DO  420 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN1).EQ.0) GOTO 420
C
      IRPC   = DIRPRD(IRPE,IRPEC)
      IRPEI = DIRPRD(IRPE,IRPI)
      IRPAB = DIRPRD(IRREPW,IRPEI)
C
      NEED = IINTFP * I030
      I020 = I010 + IRPDPD(IRPAB,ISPIN1) * VRT(IRPE,ISPIN1)
      NEED = IINTFP * I020
C
      IF(NEED.GT.MAXCOR)THEN
        WRITE(6,1010)
        CALL INSMEM('GT3WT23',NEED,MAXCOR)
      ENDIF
C
      CALL GETLST(CORE(I010),
     &       IOFFVO(IRPI,IRPEI,ISPIN1) + (I-1)*VRT(IRPE,ISPIN1) + 1,
     &            VRT(IRPE,ISPIN1),2,IRPEI,LISWVB)
C
      CALL XGEMM('N','N',
     &           IRPDPD(IRPAB,ISPIN1),VRT(IRPC,ISPIN2),VRT(IRPE,ISPIN1),
     &           ONEM,
     &           CORE(I010),IRPDPD(IRPAB,ISPIN1),
     &           CORE(I000 + IADT2(IRPC) - 1),VRT(IRPE,ISPIN1),
     &           ONE,
     &           T3(IADT3(IRPC)),IRPDPD(IRPAB,ISPIN1))
C  
  420 CONTINUE
C
C
C-----------------------------------------------------------------------
C
C     T(Cb,Km) * W(ij,ma) - T(Ca,Km) * W(ij,mb)
C
C-----------------------------------------------------------------------
C
      IRPMA = DIRPRD(IRREPW,IRPIJ)
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPMA, 8+ISPIN1)
C
      IF(IRPI.EQ.IRPJ)THEN
        IJ = IOFFOO(IRPJ,IRPIJ,ISPIN1) + INDEX(J-1) + I
      ELSE
        IJ = IOFFOO(IRPJ,IRPIJ,ISPIN1) + (J-1) * POP(IRPI,ISPIN1) + I
      ENDIF
      CALL GETLST(CORE(I000),IJ,1,2,IRPIJ,LISWOB)
C
      DO  510 IRPM=1,NIRREP
        IRPMK = DIRPRD(IRPM,IRPK)
        IRPBC = DIRPRD(IRREPT,IRPMK)
        LENT2(IRPM) = IRPDPD(IRPBC,13) * POP(IRPM,ISPIN1)
          IF(IRPM.EQ.1)THEN
            IADT2(IRPM) = 1
          ELSE
            IADT2(IRPM) = IADT2(IRPM-1) + LENT2(IRPM-1)
          ENDIF
  510 CONTINUE
C
      IOFFV = 0
      DO  520 IRPA=1,NIRREP
C
      IRPM  = DIRPRD(IRPA,IRPMA)
      IRPMK = DIRPRD(IRPM,IRPK)
      IRPBC = DIRPRD(IRREPT,IRPMK)
C
      IF(POP(IRPM,ISPIN1).EQ.0.OR.VRT(IRPA,ISPIN1).EQ.0) GOTO 520
C
      I020 = I010 + IRPDPD(IRPBC,13) * POP(IRPM,ISPIN1)
      NEED = I020 * IINTFP 
C
      IF(NEED.GT.MAXCOR)THEN
        WRITE(6,1010)
        CALL INSMEM('GT3WT23',NEED,MAXCOR)
      ENDIF
C
      DO  515 M=1,POP(IRPM,ISPIN1)
      CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPBC,13)),
     &            IOFFOO(IRPM,IRPMK,5) + (M-1)*POP(IRPK,ISPIN2) + K,
     &            1,1,IRPMK,LIST2OFF + 3)
  515 CONTINUE
C
C     T matrix is ordered (Cb,m), V matrix is (m,a).
C
      CALL XGEMM('T','T',
     &           VRT(IRPA,ISPIN1),IRPDPD(IRPBC,13),POP(IRPM,ISPIN1),
     &           ONE,
     &           CORE(I000 + IOFFV),POP(IRPM,ISPIN1),
     &           CORE(I010)        ,IRPDPD(IRPBC,13),ONE,
     &           W(IADW(IRPBC)),VRT(IRPA,ISPIN1))
C
      IOFFV = IOFFV + POP(IRPM,ISPIN1) * VRT(IRPA,ISPIN1)
  520 CONTINUE
C
C
C-----------------------------------------------------------------------
C
C     T(Cb,Mi) * W(Kj,Ma) - T(Ca,Mi) * W(Kj,Mb)
C
C-----------------------------------------------------------------------
C
C     Integrals are from list LWIC14 and are (M,a), and the
C     record is Kj.
C
      IRPMA = DIRPRD(IRREPW,IRPJK)
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPMA,12)
C
      KJ = IOFFOO(IRPJ,IRPJK,5) + (J-1)*POP(IRPK,ISPIN2) + K
      CALL GETLST(CORE(I000),KJ,1,2,IRPJK,LISWOD)
C
      I020 = I010 + IRPDPD(IRPMA,12)
      I030 = I020 + IRPDPD(IRPMA,12)
      CALL SYMTR3(IRPMA,POP(1,1),VRT(1,2),IRPDPD(IRPMA,12),1,CORE(I000),
     &            CORE(I010),CORE(I020),CORE(I030))
C
        IOFFV = 0
        DO  620 IRPM=1,NIRREP
C
        IRPMI = DIRPRD(IRPM,IRPI)
        IRPBC = DIRPRD(IRREPT,IRPMI)
        LENT2(IRPM) = IRPDPD(IRPBC,13) * POP(IRPM,ISPIN2)
          IF(IRPM.EQ.1)THEN
            IADT2(IRPM) = 1
          ELSE
            IADT2(IRPM) = IADT2(IRPM-1) + LENT2(IRPM-1)
          ENDIF
C
        IRPA = DIRPRD(IRPM,IRPMA)
        IF(POP(IRPM,ISPIN2).EQ.0.OR.VRT(IRPA,ISPIN1).EQ.0) GOTO 620
C
        I020 = I010 + IRPDPD(IRPBC,13) * POP(IRPM,ISPIN2)
        NEED = IINTFP * I020
C
          IF(NEED.GT.MAXCOR)THEN
            WRITE(6,1010)
            CALL INSMEM('GT3WT22',NEED,MAXCOR)
          ENDIF
C
        DO  610 M=1,POP(IRPM,ISPIN2)
        CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPBC,13)),
     &              IOFFOO(IRPI,IRPMI,5) + (I-1)*POP(IRPM,ISPIN2) + M,
     &              1,1,IRPMI,LIST2OFF + 3)
  610   CONTINUE
C
        CALL XGEMM('N','T',
     &             VRT(IRPA,ISPIN1),IRPDPD(IRPBC,13),POP(IRPM,ISPIN2),
     &             ONE,
     &             CORE(I000 + IOFFV),VRT(IRPA,ISPIN1),
     &             CORE(I010        ),IRPDPD(IRPBC,13),ONE,
     &             W(IADW(IRPBC)),VRT(IRPA,ISPIN1))
C
        IOFFV = IOFFV + VRT(IRPA,ISPIN1) * POP(IRPM,ISPIN2)
  620   CONTINUE
C
C-----------------------------------------------------------------------
C
C    - T(Cb,Mj) * W(Ki,Ma) + T(Ca,Mj) * W(Ki,Mb)
C
C-----------------------------------------------------------------------
C
C     Integrals are from list LWIC14 and are (M,a), and the
C     record is Ki.
C
      IRPMA = DIRPRD(IRREPW,IRPIK)
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPMA,12)
C
      KI = IOFFOO(IRPI,IRPIK,5) + (I-1)*POP(IRPK,ISPIN2) + K
      CALL GETLST(CORE(I000),KI,1,2,IRPIK,LISWOD)
C
      I020 = I010 + IRPDPD(IRPMA,12)
      I030 = I020 + IRPDPD(IRPMA,12)
      CALL SYMTR3(IRPMA,POP(1,1),VRT(1,2),IRPDPD(IRPMA,12),1,CORE(I000),
     &            CORE(I010),CORE(I020),CORE(I030))
C
        IOFFV = 0
        DO  720 IRPM=1,NIRREP
C
        IRPMJ = DIRPRD(IRPM,IRPJ)
        IRPBC = DIRPRD(IRREPT,IRPMJ)
        LENT2(IRPM) = IRPDPD(IRPBC,13) * POP(IRPM,ISPIN2)
          IF(IRPM.EQ.1)THEN
            IADT2(IRPM) = 1
          ELSE
            IADT2(IRPM) = IADT2(IRPM-1) + LENT2(IRPM-1)
          ENDIF
C
        IRPA = DIRPRD(IRPM,IRPMA)
        IF(POP(IRPM,ISPIN2).EQ.0.OR.VRT(IRPA,ISPIN1).EQ.0) GOTO 720
C
        I020 = I010 + IRPDPD(IRPBC,13) * POP(IRPM,ISPIN2)
        NEED = IINTFP * I020
C
          IF(NEED.GT.MAXCOR)THEN
            WRITE(6,1010)
            CALL INSMEM('GT3WT23',NEED,MAXCOR)
          ENDIF
C
        DO  710 M=1,POP(IRPM,ISPIN2)
        CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPBC,13)),
     &              IOFFOO(IRPJ,IRPMJ,5) + (J-1)*POP(IRPM,ISPIN2) + M,
     &              1,1,IRPMJ,LIST2OFF + 3)
  710   CONTINUE
C
        CALL XGEMM('N','T',
     &             VRT(IRPA,ISPIN1),IRPDPD(IRPBC,13),POP(IRPM,ISPIN2),
     &             ONEM,
     &             CORE(I000 + IOFFV),VRT(IRPA,ISPIN1),
     &             CORE(I010        ),IRPDPD(IRPBC,13),ONE,
     &             W(IADW(IRPBC)),VRT(IRPA,ISPIN1))
C
        IOFFV = IOFFV + VRT(IRPA,ISPIN1) * POP(IRPM,ISPIN2)
  720   CONTINUE
C
C-----------------------------------------------------------------------
C
C     -T(ab,im) * W(Kj,Cm)
C
C-----------------------------------------------------------------------
C
      IRPMC = DIRPRD(IRREPW,IRPJK)
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPMC,11)
      KJ = IOFFOO(IRPJ,IRPJK,5) + (J-1)*POP(IRPK,ISPIN2) + K
      CALL GETLST(CORE(I000),KJ,1,2,IRPJK,LISWOC)
C
      I020 = I010 + IRPDPD(IRPMC,11)
      I030 = I020 + IRPDPD(IRPMC,11)
      CALL SYMTR3(IRPMC,VRT(1,1),POP(1,2),IRPDPD(IRPMC,11),1,CORE(I000),
     &            CORE(I010),CORE(I020),CORE(I030))
C
      IOFFV = 0
      DO  850 IRPC=1,NIRREP
C
      IRPM = DIRPRD(IRPC,IRPMC)
C
      IRPIM = DIRPRD(IRPM,IRPI)
      IRPAB = DIRPRD(IRREPT,IRPIM)
C
      IF(POP(IRPM,ISPIN1).EQ.0.OR.VRT(IRPC,ISPIN2).EQ.0) GOTO 850
      IF(IRPDPD(IRPIM,ISPIN1+2).EQ.0)                    GOTO 850
C
      I020 = I010 + IRPDPD(IRPAB,ISPIN1) * POP(IRPM,ISPIN1)
      NEED = IINTFP * I020
      IF(NEED.GT.MAXCOR)THEN
        WRITE(6,1010)
        CALL INSMEM('GT3WT23',NEED,MAXCOR)
      ENDIF
C
        IF(IRPI.GT.IRPM)THEN
          IMOFF = IOFFOO(IRPI,IRPIM,ISPIN1) + (I-1)*POP(IRPM,ISPIN1) + 1
          CALL GETLST(CORE(I010),IMOFF,POP(IRPM,ISPIN1),1,IRPIM,
     &                LIST2OFF + ISPIN1)
          CALL VMINUS(CORE(I010),IRPDPD(IRPAB,ISPIN1)*POP(IRPM,ISPIN1))
        ENDIF
C
        IF(IRPI.LT.IRPM)THEN
          DO  820 M=1,POP(IRPM,ISPIN1)
          IM = IOFFOO(IRPM,IRPIM,ISPIN1) + (M-1)*POP(IRPI,ISPIN1) + I
          CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),
     &                IM,1,1,IRPIM,LIST2OFF + ISPIN1)
  820     CONTINUE
        ENDIF
C
        IF(IRPI.EQ.IRPM)THEN
          DO  830 M=1,POP(IRPM,ISPIN1)
C
            IF(I.GT.M)THEN
              IM = IOFFOO(IRPI,IRPIM,ISPIN1) + INDEX(I-1) + M
              CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),IM,
     &                    1,1,IRPIM,LIST2OFF + ISPIN1)
              CALL VMINUS(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),
     &                    IRPDPD(IRPAB,ISPIN1))
            ENDIF
C
            IF(I.LT.M)THEN
              IM = IOFFOO(IRPM,IRPIM,ISPIN1) + INDEX(M-1) + I
              CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),IM,
     &                    1,1,IRPIM,LIST2OFF + ISPIN1)
             ENDIF
C
            IF(I.EQ.M)THEN
              CALL   ZERO(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),
     &                    IRPDPD(IRPAB,ISPIN1))
            ENDIF
  830     CONTINUE
        ENDIF
C
      CALL XGEMM('N','N',
     &           IRPDPD(IRPAB,ISPIN1),VRT(IRPC,ISPIN2),POP(IRPM,ISPIN1),
     &           ONEM,
     &           CORE(I010),IRPDPD(IRPAB,ISPIN1),
     &           CORE(I000 + IOFFV),POP(IRPM,ISPIN1),
     &           ONE,
     &           T3(IADT3(IRPC)),IRPDPD(IRPAB,ISPIN1))
C
      IOFFV = IOFFV + POP(IRPM,ISPIN1) * VRT(IRPC,ISPIN2)
  850 CONTINUE
C
C
C-----------------------------------------------------------------------
C
C      T(ab,jm) * W(Ki,Cm)
C
C-----------------------------------------------------------------------
C
      IRPMC = DIRPRD(IRREPW,IRPIK)
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPMC,11)
      KI = IOFFOO(IRPI,IRPIK,5) + (I-1)*POP(IRPK,ISPIN2) + K
      CALL GETLST(CORE(I000),KI,1,2,IRPIK,LISWOC)
C
      I020 = I010 + IRPDPD(IRPMC,11)
      I030 = I020 + IRPDPD(IRPMC,11)
      CALL SYMTR3(IRPMC,VRT(1,1),POP(1,2),IRPDPD(IRPMC,11),1,CORE(I000),
     &            CORE(I010),CORE(I020),CORE(I030))
C
      IOFFV = 0
      DO  950 IRPC=1,NIRREP
C
      IRPM = DIRPRD(IRPC,IRPMC)
C
      IRPJM = DIRPRD(IRPM,IRPJ)
      IRPAB = DIRPRD(IRREPT,IRPJM)
C
      IF(POP(IRPM,ISPIN1).EQ.0.OR.VRT(IRPC,ISPIN2).EQ.0) GOTO 950
C
      I020 = I010 + IRPDPD(IRPAB,ISPIN1) * POP(IRPM,ISPIN1)
      NEED = IINTFP * I020
      IF(NEED.GT.MAXCOR)THEN
        WRITE(6,1010)
        CALL INSMEM('GT3WT23',NEED,MAXCOR)
      ENDIF
C
      IF(IRPJ.GT.IRPM)THEN
        JMOFF = IOFFOO(IRPJ,IRPJM,ISPIN1) + (J-1)*POP(IRPM,ISPIN1) + 1
        CALL GETLST(CORE(I010),JMOFF,POP(IRPM,ISPIN1),1,IRPJM,
     &              LIST2OFF + ISPIN1)
        CALL VMINUS(CORE(I010),IRPDPD(IRPAB,ISPIN1)*POP(IRPM,ISPIN1))
      ENDIF
C
      IF(IRPJ.LT.IRPM)THEN
        DO  920 M=1,POP(IRPM,ISPIN1)
        JM = IOFFOO(IRPM,IRPJM,ISPIN1) + (M-1)*POP(IRPJ,ISPIN1) + J
        CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),
     &              JM,1,1,IRPJM,LIST2OFF + ISPIN1)
  920   CONTINUE
      ENDIF
C
      IF(IRPJ.EQ.IRPM)THEN
        DO  930 M=1,POP(IRPM,ISPIN1)
C
          IF(J.GT.M)THEN
            JM = IOFFOO(IRPJ,IRPJM,ISPIN1) + INDEX(J-1) + M
            CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),JM,
     &                  1,1,IRPJM,LIST2OFF + ISPIN1)
            CALL VMINUS(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),
     &                  IRPDPD(IRPAB,ISPIN1))
          ENDIF
C
          IF(J.LT.M)THEN
            JM = IOFFOO(IRPM,IRPJM,ISPIN1) + INDEX(M-1) + J
            CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),JM,
     &                  1,1,IRPJM,LIST2OFF + ISPIN1)
          ENDIF
C
          IF(J.EQ.M)THEN
            CALL   ZERO(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),
     &                  IRPDPD(IRPAB,ISPIN1))
          ENDIF
  930   CONTINUE
        ENDIF
C
      CALL XGEMM('N','N',
     &           IRPDPD(IRPAB,ISPIN1),VRT(IRPC,ISPIN2),POP(IRPM,ISPIN1),
     &           ONE,
     &           CORE(I010),IRPDPD(IRPAB,ISPIN1),
     &           CORE(I000 + IOFFV),POP(IRPM,ISPIN1),
     &           ONE,
     &           T3(IADT3(IRPC)),IRPDPD(IRPAB,ISPIN1))
C
      IOFFV = IOFFV + POP(IRPM,ISPIN1) * VRT(IRPC,ISPIN2)
  950 CONTINUE
C
      RETURN
 1010 FORMAT(' @GT3WT23-I, Insufficient memory to continue. ')
      END
