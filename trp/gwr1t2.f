      SUBROUTINE GWR1T2(CORE,MAXCOR,IUHF,IRREPT,LISTT1,
     &                  LABCI0,LLCJK0,DOABCI,DOLCJK,ADDFOV,
     &                  FACABCI,FACLCJK)
C
C     This subroutine first computes the modified O-V quantity
C
C     F(f,l) = \sum_{n,g} T(g,n) * <nl||gf>
C
C     The integral is assumed to be totally symmetric, but the amplitude
C     need not be. Once F is computed, it is contracted with T2 (assumed
C     to be totally symmetric) to give contributions to W(bcdk) and W(lcjk).
C     Note that F is never stored on disk.
C
C     DOABCI ---- Whether to calculate the contribution to W(abci).
C     DOLCJK ---- Whether to calculate the contribution to W(lcjk).
C     ADDFOV ---- Whether to add the O-V block of the Fock matrix
C                 in the non-HF case.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT,DSZW,D,DK
      LOGICAL DOABCI,DOLCJK,ADDFOV
      DIMENSION CORE(1),I0F(2),I0T1(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFEA(2),NFMI(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON/FLAGS/ IFLAGS(100)
C
      DATA ONE,ONEM,ZILCH,TWO,HALF,HALFM 
     &     /1.0D0,-1.0D0,0.0D0,2.0D0,0.5D0,-0.5D0/
C
      INDEX(I) = I*(I-1) / 2
C
      I0F(1)  = 1
      I0F(2)  = I0F(1)  + IRPDPD(IRREPT, 9)
      I0T1(1) = I0F(2)  + IRPDPD(IRREPT,10)
      I0T1(2) = I0T1(1) + IRPDPD(IRREPT, 9)
      I000    = I0T1(2) + IRPDPD(IRREPT,10)
C
      CALL ZERO(CORE(I0F(1)) ,IRPDPD(IRREPT, 9))
      CALL ZERO(CORE(I0F(2)) ,IRPDPD(IRREPT,10))
      CALL ZERO(CORE(I0T1(1)),IRPDPD(IRREPT, 9))
      CALL ZERO(CORE(I0T1(2)),IRPDPD(IRREPT,10))
C
C     ----- Get the T1 vectors -----
C
      CALL GETLST(CORE(I0T1(1)),1,1,1,1     ,LISTT1)
      CALL GETLST(CORE(I0T1(2)),1,1,1,1+IUHF,LISTT1)
C
      IRPGN = IRREPT
      IRPFL = IRPGN
C
C     F(F,L) = T(G,N) * <GF || NL>
C     F(f,l) = T(g,n) * <gf || nl>
C
      DO 10 ISPIN=1,1+IUHF
C
      LISTW = 18 + ISPIN
C
C     ----- Integral is stored (GN;FL) -----
C
      DSZW = IRPDPD(IRPGN,8+ISPIN)
      NDSW = IRPDPD(IRPFL,8+ISPIN)
      I010 = I000 + DSZW * NDSW
C
      CALL GETLST(CORE(I000),1,NDSW,2,IRPFL,LISTW)
C
C     F(FL) = (GN;FL)^T * T1(GN)
C
      CALL XGEMM('T','N',NDSW,1,DSZW,ONE,CORE(I000),DSZW,
     &           CORE(I0T1(ISPIN)),DSZW,ONE,CORE(I0F(ISPIN)),NDSW)
   10 CONTINUE
C
C    F(F,L) = T1(g,n) * <Ln || Fg>    (FL;gn)
C    F(f,l) = T1(G,N) * <Nl || Gf>    (GN;fl)
C 
      IRPGN = IRREPT
      IRPFL = IRPGN
C
      DSZW = IRPDPD(IRPFL, 9)
      NDSW = IRPDPD(IRPGN,10)
C
      LISTW = 18
C
      CALL GETLST(CORE(I000),1,NDSW,2,IRPGN,LISTW)
C
C     F(FL) = (FL;gn) * T1(gn)
C
      CALL XGEMM('N','N',DSZW,1,NDSW,ONE,CORE(I000),DSZW,
     &           CORE(I0T1(2)),NDSW,ONE,CORE(I0F(1)),DSZW)
C
      IF(IUHF.NE.0)THEN
C
C     F(fl) = (GN;fl)^T * T1(GN)
C
      CALL XGEMM('T','N',NDSW,1,DSZW,ONE,CORE(I000),DSZW,
     &           CORE(I0T1(1)),DSZW,ONE,CORE(I0F(2)),NDSW)
      ENDIF
C
C     Add Fock matrix if this has been requested. Current code assumes
C     the Fock matrix is just that, ie not some derivative, hence specific
C     list number.
C
      IF(ADDFOV)THEN
        IF(IFLAGS(38).GT.0)THEN
          DO 15 ISPIN=1,IUHF+1
            CALL GETLST(CORE(I000),1,1,2,2+ISPIN,93)
            CALL  SAXPY(IRPDPD(IRREPT,8+ISPIN),ONE,CORE(I000),1,
     &                  CORE(I0F(ISPIN)),1)
   15     CONTINUE
        ENDIF
      ENDIF
C
      IF(IUHF.EQ.0)THEN
c YAU : old
c       CALL ICOPY(DSZW*IINTFP,CORE(I0F(1)),1,CORE(I0F(2)),1)
c YAU : new
        CALL DCOPY(DSZW,CORE(I0F(1)),1,CORE(I0F(2)),1)
c YAU : end
      ENDIF
C
C     ----- Now that we have the F intermediates, we can compute the -----
C     ----- contributions to W(lcjk) and W(bcdk)                     -----
C
      IF(DOLCJK)THEN
C
C     W(lcjk) = - \sum_f   T(cf,jk) F(f,l)
C
C         LCJK        CFJK    FL
C         lcjk        cfjk    fl
C
      IF(IUHF.NE.0)THEN
C
      DO 50 ISPIN=1,2
C
      DO 40 IRPLC=1,NIRREP
      IRPJK = DIRPRD(IRPLC,IRREPT)
C
      DSZW = IRPDPD(IRPJK,2 + ISPIN)
      NDSW = IRPDPD(IRPLC,8 + ISPIN)
C
      I010 = I000 + DSZW * NDSW
      I020 = I010 + DSZW * NDSW
C
      CALL GETLST(CORE(I010),1,NDSW,2,IRPLC,LLCJK0+ISPIN)
C
C     Reorder this symmetry of W from (JK,LC) to (LC,JK).
C
      CALL TRANSP(CORE(I010),CORE(I000),NDSW,DSZW)
C
      DO 30 JK=1,DSZW
C
      CALL GETLST(CORE(I010),JK,1,1,IRPJK,43 + ISPIN)
      IRPCF = IRPJK
      CALL SYMEXP2(IRPCF,VRT(1,ISPIN),IRPDPD(IRPCF,18+ISPIN),
     &             IRPDPD(IRPCF,ISPIN),IRPDPD(IRPJK,2+ISPIN),
     &             CORE(I010),CORE(I010))
C
      DO 20 IRPF=1,NIRREP
C
      IRPC  = DIRPRD(IRPF,IRPCF)
      IRPL  = DIRPRD(IRREPT,IRPF)
C
      IOFFT = I010 + ISYMOFF(IRPF,IRPCF,18+ISPIN) - 1
      IOFFF = I0F(ISPIN) + ISYMOFF(IRPL,IRREPT,8+ISPIN) - 1
      IOFFW = I000 + (JK-1)*NDSW + ISYMOFF(IRPC,IRPLC,15+ISPIN) - 1
C
C     W(L,C) = W(L,C) - F(F,L) * T(C,F)
C
      CALL XGEMM('T','T',
     &           POP(IRPL,ISPIN),VRT(IRPC,ISPIN),VRT(IRPF,ISPIN),
C    &           HALFM,
     &           -FACLCJK,
     &           CORE(IOFFF),VRT(IRPF,ISPIN),
     &           CORE(IOFFT),VRT(IRPC,ISPIN),ONE,
     &           CORE(IOFFW),POP(IRPL,ISPIN))
   20 CONTINUE
   30 CONTINUE
C
C     (LC,JK) ---> (JK,LC).
C
      CALL TRANSP(CORE(I000),CORE(I010),DSZW,NDSW)
      CALL PUTLST(CORE(I010),1,NDSW,2,IRPLC,LLCJK0+ISPIN)
   40 CONTINUE
   50 CONTINUE
C
C         lCJk        CfJk    fl    (W storage is Cl)
C
      DO 90 IRPLC=1,NIRREP
      IRPCL = IRPLC
      IRPJK = DIRPRD(IRPLC,IRREPT)
C
      DSZW = IRPDPD(IRPJK,14)
      NDSW = IRPDPD(IRPLC,11)
C
      I010 = I000 + DSZW * NDSW
      I020 = I010 + DSZW * NDSW
C
      CALL GETLST(CORE(I010),1,NDSW,2,IRPLC,LLCJK0+3)
C
C     Reorder this symmetry of W from (Jk,Cl) to (Cl,Jk).
C
      CALL TRANSP(CORE(I010),CORE(I000),NDSW,DSZW)
C
      DO 80 JK=1,DSZW
C
      CALL GETLST(CORE(I010),JK,1,1,IRPJK,46)
      IRPCF = IRPJK
C
      DO 70 IRPF=1,NIRREP
C
      IRPC  = DIRPRD(IRPF,IRPCF)
      IRPL  = DIRPRD(IRREPT,IRPF)
C
      IOFFT = I010   +               ISYMOFF(IRPF,IRPCF ,13) - 1
      IOFFF = I0F(2) +               ISYMOFF(IRPL,IRREPT,10) - 1
      IOFFW = I000   + (JK-1)*NDSW + ISYMOFF(IRPL,IRPLC ,11) - 1
C
C     W(C,l) = W(C,l) - T(C,f) * F(f,l). This gives spin case JklC. To
C     get JkCl we negate, i.e. put a plus in the above.
C
      CALL XGEMM('N','N',
     &           VRT(IRPC,1),POP(IRPL,2),VRT(IRPF,2),
C    &           HALF,
     &           FACLCJK,
     &           CORE(IOFFT),VRT(IRPC,1),
     &           CORE(IOFFF),VRT(IRPF,2),ONE,
     &           CORE(IOFFW),VRT(IRPC,1))
   70 CONTINUE
   80 CONTINUE
C
C     (Cl,Jk) ---> (Jk,Cl).
C
      CALL TRANSP(CORE(I000),CORE(I010),DSZW,NDSW)
      CALL PUTLST(CORE(I010),1,NDSW,2,IRPLC,LLCJK0+3)
   90 CONTINUE
C
      ENDIF
C
C         LcJk        cFJk    FL
C        W(LcJk) = -T(cF,Jk) * F(F,L) = T(Fc,Jk) * F(F,L)
C
      DO 150 IRPLC=1,NIRREP
      IRPJK = DIRPRD(IRPLC,IRREPT)
C
      DSZW = IRPDPD(IRPJK,14)
      NDSW = IRPDPD(IRPLC,18)
C
      IF(DSZW.EQ.0.OR.NDSW.EQ.0) GOTO 150
C
      I010 = I000 + DSZW * NDSW
      I020 = I010 + DSZW * NDSW
C
      CALL GETLST(CORE(I010),1,NDSW,2,IRPLC,LLCJK0+4)
C
C     Reorder this symmetry of W from (Jk,Lc) to (Lc,Jk).
C
      CALL TRANSP(CORE(I010),CORE(I000),NDSW,DSZW)
C
      DO 140 JK=1,DSZW
C
      CALL GETLST(CORE(I010),JK,1,1,IRPJK,46)
      IRPFC = IRPJK
C
      DO 130 IRPF=1,NIRREP
C
      IRPC  = DIRPRD(IRPF,IRPFC)
      IRPL  = DIRPRD(IRREPT,IRPF)
C
      IF(POP(IRPL,1).EQ.0.OR.VRT(IRPC,2).EQ.0.OR.VRT(IRPF,1).EQ.0)
     &   GOTO 130
C
      IOFFT = I010   +               ISYMOFF(IRPC,IRPFC ,13) - 1
      IOFFF = I0F(1) +               ISYMOFF(IRPL,IRREPT, 9) - 1
      IOFFW = I000   + (JK-1)*NDSW + ISYMOFF(IRPC,IRPLC ,18) - 1
C
C     W(L,c) = W(L,c) + F(F,L) * T(F,c)
C
      CALL XGEMM('T','N',
     &           POP(IRPL,1),VRT(IRPC,2),VRT(IRPF,1),
C    &           HALF,
     &           FACLCJK,
     &           CORE(IOFFF),VRT(IRPF,1),
     &           CORE(IOFFT),VRT(IRPF,1),ONE,
     &           CORE(IOFFW),POP(IRPL,1))
  130 CONTINUE
  140 CONTINUE
C
C     (Lc,Jk) ---> (Jk,Lc).
C
      CALL TRANSP(CORE(I000),CORE(I010),DSZW,NDSW)
      CALL PUTLST(CORE(I010),1,NDSW,2,IRPLC,LLCJK0+4)
  150 CONTINUE
C
      ENDIF
C
      IF(DOABCI)THEN
C
C     W(bcdk) = \sum_n F(d,n) * T(bc,kn)
C
C       BcDk    DN    BckN. Evaluate as -F(D,N)*T(Bc,Nk).
C
      DO 250 IRPDK=1,NIRREP
      IRPBC = DIRPRD(IRREPT,IRPDK)
      IRPNK = IRPBC
C
      NDK = IRPDPD(IRPDK,11)
      NBC = IRPDPD(IRPBC,13)
      IF(NDK.EQ.0.OR.NBC.EQ.0) GOTO 250
C
      DO 240 IRPK =1,NIRREP
      IRPD  = DIRPRD(IRPK,IRPDK)
      IRPN  = DIRPRD(IRREPT,IRPD)
C
      IF(POP(IRPK,2).EQ.0.OR.VRT(IRPD,1).EQ.0.OR.POP(IRPN,1).EQ.0)
     &                         GOTO 240
C
      DO 230 K=1,POP(IRPK,2)
      DO 220 D=1,VRT(IRPD,1)
C
      DK = ISYMOFF(IRPK,IRPDK,11) - 1 + (K-1)*VRT(IRPD,1) + D
      CALL GETLST(CORE(I000),DK,1,2,IRPDK,LABCI0+4)
C
      I010 = I000 + IRPDPD(IRPBC,13)
C
      DO 210 N=1,POP(IRPN,1)
      NK = ISYMOFF(IRPK,IRPNK,14) - 1 + (K-1)*POP(IRPN,1) + N
      CALL GETLST(CORE(I010),NK,1,1,IRPNK,46)
C
      IOFFF = I0F(1) + ISYMOFF(IRPN,IRREPT,9) - 1 + 
     &                           (N-1)*VRT(IRPD,1) + D - 1
C
C     FACT = HALFM * CORE(IOFFF)
      FACT = -FACABCI * CORE(IOFFF)
      CALL SAXPY(IRPDPD(IRPBC,13),FACT,CORE(I010),1,CORE(I000),1)
  210 CONTINUE
      CALL PUTLST(CORE(I000),DK,1,2,IRPDK,LABCI0+4)
  220 CONTINUE
  230 CONTINUE
  240 CONTINUE
  250 CONTINUE
C
        IF(IUHF.NE.0)THEN
C
C     W(bcdk) =   \sum_n F(d,n) * T(bc,kn)
C     W(bckd) = - \sum_n F(d,n) * T(bc,kn)
C
C       BcKd    dn    BcKn. Evaluate as -F(d,n)*T(Bc,Kn).
C
          DO 350 IRPKD=1,NIRREP
          IRPBC = DIRPRD(IRREPT,IRPKD)
          IRPKN = IRPBC
C
          NKD = IRPDPD(IRPKD,18)
          NBC = IRPDPD(IRPBC,13)
          IF(NKD.EQ.0.OR.NBC.EQ.0) GOTO 350
C
          DO 340 IRPK =1,NIRREP
          IRPD  = DIRPRD(IRPK,IRPKD)
          IRPN  = DIRPRD(IRREPT,IRPD)
C
          IF(POP(IRPK,1).EQ.0.OR.VRT(IRPD,2).EQ.0.OR.POP(IRPN,2).EQ.0)
     &                             GOTO 340
C
          DO 330 K=1,POP(IRPK,1)
          DO 320 D=1,VRT(IRPD,2)
C
          KD = ISYMOFF(IRPD,IRPKD,18) - 1 + (D-1)*POP(IRPK,1) + K
          CALL GETLST(CORE(I000),KD,1,2,IRPKD,LABCI0+3)
C
          I010 = I000 + IRPDPD(IRPBC,13)
C
          DO 310 N=1,POP(IRPN,1)
          KN = ISYMOFF(IRPN,IRPKN,14) - 1 + (N-1)*POP(IRPK,1) + K
          CALL GETLST(CORE(I010),KN,1,1,IRPKN,46)
C
          IOFFF = I0F(2) + ISYMOFF(IRPN,IRREPT,10) - 1 + 
     &                               (N-1)*VRT(IRPD,2) + D - 1
C
C         FACT = HALFM * CORE(IOFFF)
          FACT = -FACABCI * CORE(IOFFF)
          CALL SAXPY(IRPDPD(IRPBC,13),FACT,CORE(I010),1,CORE(I000),1)
  310     CONTINUE
          CALL PUTLST(CORE(I000),KD,1,2,IRPKD,LABCI0+3)
  320     CONTINUE
  330     CONTINUE
  340     CONTINUE
  350     CONTINUE
C
C
C     W(bcdk) =   \sum_n F(d,n) * T(bc,kn)
C
C       BCDK    DN    BCKN. Evaluate as +F(D,N)*T(BC,KN).
C
          DO 460 ISPIN=1,2
          DO 450 IRPDK=1,NIRREP
          IRPBC = DIRPRD(IRREPT,IRPDK)
          IRPKN = IRPBC
C
          NDK = IRPDPD(IRPDK,8+ISPIN)
          NBC = IRPDPD(IRPBC,ISPIN)
          IF(NDK.EQ.0.OR.NBC.EQ.0) GOTO 450
C
          DO 440 IRPK =1,NIRREP
          IRPD  = DIRPRD(IRPK,IRPDK)
          IRPN  = DIRPRD(IRREPT,IRPD)
C
          IF(POP(IRPK,ISPIN).EQ.0.OR.VRT(IRPD,ISPIN).EQ.0
     &                           .OR.POP(IRPN,ISPIN).EQ.0) GOTO 440
C
          DO 430 K=1,POP(IRPK,ISPIN)
          DO 420 D=1,VRT(IRPD,ISPIN)
C
          DK = ISYMOFF(IRPK,IRPDK,8+ISPIN) - 1 + (K-1)*VRT(IRPD,ISPIN)
     &                                          + D
          CALL GETLST(CORE(I000),DK,1,2,IRPDK,LABCI0+ISPIN)
C
          I010 = I000 + IRPDPD(IRPBC,ISPIN)
C
          DO 410 N=1,POP(IRPN,ISPIN)
C
          IF(IRPK.LT.IRPN)THEN
            KN = ISYMOFF(IRPN,IRPKN,2+ISPIN) - 1 + (N-1)*POP(IRPK,ISPIN)
     &                                           +  K
            CALL GETLST(CORE(I010),KN,1,1,IRPKN,43+ISPIN)
          ELSEIF(IRPK.GT.IRPN)THEN
            KN = ISYMOFF(IRPK,IRPKN,2+ISPIN) - 1 + (K-1)*POP(IRPN,ISPIN)
     &                                           +  N
            CALL GETLST(CORE(I010),KN,1,1,IRPKN,43+ISPIN)
            CALL VMINUS(CORE(I010),NBC)
          ELSEIF(IRPK.EQ.IRPN)THEN
            IF(K.LT.N)THEN
              KN = ISYMOFF(IRPN,IRPKN,2+ISPIN) - 1 + INDEX(N-1) + K
              CALL GETLST(CORE(I010),KN,1,1,IRPKN,43+ISPIN)
            ELSEIF(K.GT.N)THEN
              KN = ISYMOFF(IRPK,IRPKN,2+ISPIN) - 1 + INDEX(K-1) + N
              CALL GETLST(CORE(I010),KN,1,1,IRPKN,43+ISPIN)
              CALL VMINUS(CORE(I010),NBC)
            ELSEIF(K.EQ.N)THEN
              CALL ZERO(CORE(I010),NBC)
            ENDIF
          ENDIF
C
          IOFFF = I0F(ISPIN) + ISYMOFF(IRPN,IRREPT,8+ISPIN) - 1 + 
     &                               (N-1)*VRT(IRPD,ISPIN) + D - 1
C
          FACT =  FACABCI * CORE(IOFFF)
          CALL SAXPY(IRPDPD(IRPBC,ISPIN),FACT,CORE(I010),1,CORE(I000),1)
  410     CONTINUE
          CALL PUTLST(CORE(I000),DK,1,2,IRPDK,LABCI0+ISPIN)
  420     CONTINUE
  430     CONTINUE
  440     CONTINUE
  450     CONTINUE
  460     CONTINUE
C
        ENDIF
C
      ENDIF
C
      RETURN
      END
