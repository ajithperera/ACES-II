      SUBROUTINE GCC3W4TAU(CORE,MAXCOR,IUHF,IRREPX,
     &                     LZOFF,LWOFF,LR2OFF,LT1,LT1OFF,LR1,LR1OFF,
     &                     INCR2)
C-----------------------------------------------------------------------
C
C This subroutine is used to compute the contribution
C
C     Z(jk,lc) = - \sum_ef [ tau(r,t) (ef,jk) ] * W(clef)
C     Z        =            T * W
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      DOUBLE PRECISION CORE
      INTEGER MAXCOR,IUHF,IRREPX,LZOFF,LWOFF,LR2OFF,LT1,LT1OFF,LR1,
     &        LR1OFF
      LOGICAL INCR2
C-----------------------------------------------------------------------
      INTEGER IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      INTEGER POP,VRT,NT,NFEA,NFMI,NSTART,NIRREP,IRREPY,DIRPRD
      INTEGER IRPDPD,ISYTYP,ID,ISYMOFF
C-----------------------------------------------------------------------
      DOUBLE PRECISION ONE,ONEM,ZILCH
      INTEGER I0T1,I0R1,I000,I010,I020,I030,I040,I050,LISTW,LISTZ
      INTEGER DISSIZT,NDIST,NDISZ,C,CL,L,LC,ISPIN
      INTEGER IRPC,IRPCL,IRPEF,IRPJK,IRPL,IRPLC
C-----------------------------------------------------------------------
      DIMENSION CORE(1),I0T1(2),I0R1(2)
C-----------------------------------------------------------------------
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFEA(2),NFMI(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
C-----------------------------------------------------------------------
      DATA ONE,ONEM,ZILCH /1.0D+00,-1.0D+00,0.0D+00/
C-----------------------------------------------------------------------
C
C --- Read T1 amplitudes ---
C
      I0T1(1) = 1
      I0T1(2) = I0T1(1) + IRPDPD(     1, 9)
      CALL GETLST(CORE(I0T1(1)),1,1,1,LT1OFF + 1     ,LT1)
      CALL GETLST(CORE(I0T1(2)),1,1,1,LT1OFF + 1+IUHF,LT1)
C
C --- Read R1 amplitudes ---
C
      I0R1(1) = I0T1(2) + IRPDPD(     1,10)
      I0R1(2) = I0R1(1) + IRPDPD(IRREPX, 9)
      CALL GETLST(CORE(I0R1(1)),1,1,1,LR1OFF + 1     ,LR1)
      CALL GETLST(CORE(I0R1(2)),1,1,1,LR1OFF + 1+IUHF,LR1)
C
      I000 = I0R1(2) + IRPDPD(IRREPX,10)
C
C-----------------------------------------------------------------------
C     Z(Jk,Lc).
C-----------------------------------------------------------------------
C
C     Z(Jk,Lc) = - \sum_ef [ tau(r,t) (Ef,Jk) ] * W(cLEf)
C              =   \sum_ef [ tau(r,t) (Ef,Jk) ] * W(LcEf)
C
C     RHF :
C
C     Z(Jk,Lc) =   \sum_ef [ tau(r,t) (Ef,Jk) ] * W(cLfE)
C     W(cLfE)  = W(ClFe).
C
      LISTZ  = LZOFF + 4
      LISTW  = LWOFF + 4 - IUHF
C
      DO   100 IRPLC=1,NIRREP
C
      IRPCL = IRPLC
      IRPEF = IRPCL
      IRPJK = DIRPRD(IRREPX,IRPEF)
C
      DISSIZT = IRPDPD(IRPEF,13)
      NDIST   = IRPDPD(IRPJK,14)
      NDISZ   = IRPDPD(IRPLC,18)
C
      IF(NDIST.EQ.0 .OR. NDISZ.EQ.0 .OR. DISSIZT.EQ.0) GOTO 100
C
      I010 = I000 + DISSIZT * NDIST
      I020 = I010 + DISSIZT
      I030 = I020 + MAX(DISSIZT,NDIST)
      I040 = I030 + MAX(DISSIZT,NDIST)
      I050 = I040 + MAX(DISSIZT,NDIST)
C
      IF(INCR2)THEN
        CALL GETLST(CORE(I000),1,NDIST,1,IRPJK,LR2OFF+3)
      ELSE
        CALL ZERO(CORE(I000),DISSIZT*NDIST)
      ENDIF
      CALL DTAU(IRPEF,IRPJK,     1,IRREPX,CORE(I000),
     &          CORE(I0T1(1)),CORE(I0T1(2)),
     &          CORE(I0R1(1)),CORE(I0R1(2)),3,1.0D+00)
C
      DO    90  IRPC=1,NIRREP
C
      IRPL = DIRPRD(IRPLC,IRPC)
C
      DO  80     C=1,VRT(IRPC,2)
      DO  70     L=1,POP(IRPL,1)
C
      LC = ISYMOFF(IRPC,IRPLC,18) + (C-1)*POP(IRPL,1) + L - 1
C
      IF(IUHF.EQ.0)THEN
        CL = ISYMOFF(IRPL,IRPCL,11) + (L-1)*VRT(IRPC,2) + C - 1
        CALL GETLST(CORE(I010),CL,1,2,IRPCL,LISTW)
        CALL SYMTR3(IRPEF,VRT(1,1),VRT(1,2),IRPDPD(IRPEF,13),1,
     &              CORE(I010),CORE(I020),CORE(I030),CORE(I040))
      ELSE
        CALL GETLST(CORE(I010),LC,1,2,IRPLC,LISTW)
      ENDIF
C
      CALL GETLST(CORE(I020),LC,1,2,IRPLC,LISTZ)
C
C     Z(Jk) = TAU(Ef,Jk) * W(Ef)
      CALL XGEMM('T','N',NDIST,1,DISSIZT,
     &           ONE,
     &           CORE(I000),DISSIZT,
     &           CORE(I010),DISSIZT,ONE,
     &           CORE(I020),NDIST)
C
      CALL PUTLST(CORE(I020),LC,1,2,IRPLC,LISTZ)
   70 CONTINUE
   80 CONTINUE
   90 CONTINUE
  100 CONTINUE
C
      IF(IUHF.EQ.0) RETURN
C
C-----------------------------------------------------------------------
C     Z(Jk,Cl).
C-----------------------------------------------------------------------
C
C     Z(jk,lc) = - \sum_ef [ tau(r,t) (ef,jk) ] * W(clef)
C  => Z(Jk,Cl) =   \sum_ef [ tau(r,t) (Ef,Jk) ] * W(ClEf)
C
      LISTZ  = LZOFF + 3
      LISTW  = LWOFF + 4
C
      DO   200 IRPCL=1,NIRREP
C
      IRPEF = IRPCL
      IRPJK = DIRPRD(IRREPX,IRPEF)
C
      DISSIZT = IRPDPD(IRPEF,13)
      NDIST   = IRPDPD(IRPJK,14)
C
      NDISZ   = IRPDPD(IRPCL,11)
C
      IF(NDIST.EQ.0 .OR. NDISZ.EQ.0 .OR. DISSIZT.EQ.0) GOTO 200
C
      I010 = I000 + DISSIZT * NDIST
      I020 = I010 + DISSIZT
      I030 = I020 + MAX(DISSIZT,NDIST)
      I040 = I030 + MAX(DISSIZT,NDIST)
      I050 = I040 + MAX(DISSIZT,NDIST)
C
      IF(INCR2)THEN
        CALL GETLST(CORE(I000),1,NDIST,1,IRPJK,LR2OFF+3)
      ELSE
        CALL ZERO(CORE(I000),DISSIZT*NDIST)
      ENDIF
      CALL DTAU(IRPEF,IRPJK,     1,IRREPX,CORE(I000),
     &          CORE(I0T1(1)),CORE(I0T1(2)),
     &          CORE(I0R1(1)),CORE(I0R1(2)),3,1.0D+00)
C
      DO   190  IRPL=1,NIRREP
C
      IRPC = DIRPRD(IRPCL,IRPL)
C
      DO  180     L=1,POP(IRPL,2)
      DO  170     C=1,VRT(IRPC,1)
C
      CL = ISYMOFF(IRPL,IRPCL,11) + (L-1)*VRT(IRPC,1) + C - 1
C
      CALL GETLST(CORE(I010),CL,1,2,IRPCL,LISTW)
      CALL GETLST(CORE(I020),CL,1,2,IRPCL,LISTZ)
C
C     Z(Jk) = TAU(Ef,Jk) * W(Ef)
      CALL XGEMM('T','N',NDIST,1,DISSIZT,
     &           ONE,
     &           CORE(I000),DISSIZT,
     &           CORE(I010),DISSIZT,ONE,
     &           CORE(I020),NDIST)
C
      CALL PUTLST(CORE(I020),CL,1,2,IRPCL,LISTZ)
  170 CONTINUE
  180 CONTINUE
  190 CONTINUE
  200 CONTINUE
C
C-----------------------------------------------------------------------
C     Z(JK,LC)/Z(jk,lc).
C-----------------------------------------------------------------------
C
C     Z(jk,lc) = - \sum_ef [ tau(r,t) (ef,jk) ] * W(clef)
C
      DO   400 ISPIN=1,2
C
      LISTZ  = LZOFF + ISPIN
      LISTW  = LWOFF + ISPIN
C
      DO   300 IRPLC=1,NIRREP
C
      IRPCL = IRPLC
      IRPEF = IRPLC
      IRPJK = DIRPRD(IRREPX,IRPEF)
C
      DISSIZT = IRPDPD(IRPEF,   ISPIN)
      NDIST   = IRPDPD(IRPJK, 2+ISPIN)
      NDISZ   = IRPDPD(IRPLC,15+ISPIN)
C
      IF(DISSIZT.EQ.0 .OR. NDIST.EQ.0 .OR. NDISZ.EQ.0) GOTO 300
C
      I010 = I000 + DISSIZT * NDIST
      I020 = I010 + DISSIZT
      I030 = I020 + MAX(DISSIZT,NDIST)
      I040 = I030 + MAX(DISSIZT,NDIST)
      I050 = I040 + MAX(DISSIZT,NDIST)
C
      IF(INCR2)THEN
        CALL GETLST(CORE(I000),1,NDIST,1,IRPJK,LR2OFF+ISPIN)
      ELSE
        CALL ZERO(CORE(I000),DISSIZT*NDIST)
      ENDIF
      CALL DTAU(IRPEF,IRPJK,     1,IRREPX,CORE(I000),
     &          CORE(I0T1(1)),CORE(I0T1(2)),
     &          CORE(I0R1(1)),CORE(I0R1(2)),ISPIN,1.0D+00)
C
      DO   290  IRPC=1,NIRREP
C
      IRPL = DIRPRD(IRPLC,IRPC)
C
      DO   280     C=1,VRT(IRPC,ISPIN)
      DO   270     L=1,POP(IRPL,ISPIN)
C
      LC = ISYMOFF(IRPC,IRPLC,15+ISPIN) + (C-1)*POP(IRPL,ISPIN) + L - 1
      CL = ISYMOFF(IRPL,IRPCL, 8+ISPIN) + (L-1)*VRT(IRPC,ISPIN) + C - 1
C
      CALL GETLST(CORE(I010),CL,1,2,IRPCL,LISTW)
      CALL GETLST(CORE(I020),LC,1,2,IRPLC,LISTZ)
C
C     Z(jk) = TAU(ef,jk) * W(ef)
      CALL XGEMM('T','N',NDIST,1,DISSIZT,
     &           ONEM,
     &           CORE(I000),DISSIZT,
     &           CORE(I010),DISSIZT,ONE,
     &           CORE(I020),NDIST)
C
      CALL PUTLST(CORE(I020),LC,1,2,IRPLC,LISTZ)
  270 CONTINUE
  280 CONTINUE
  290 CONTINUE
  300 CONTINUE
  400 CONTINUE
      RETURN
      END
