      SUBROUTINE GCC3W4RING(CORE,MAXCOR,IUHF,IRREPT,
     &                      LZOFF,LWOFF,LT1,LT1OFF)
C-----------------------------------------------------------------------
C This subroutine is used to compute the contribution
C
C     Z(lc,jk) =   \sum_f [ t(f,j) * W(lcfk)  - t(f,k) * W(lcfj)  ]
C     Z        =               T   * W
C
C to the Hbar element. It is called by EOMNT3 and HBARR3. Very similar
C contractions are performed by DHBIAJK3 and CC3W4RING. CC3W4RING is spec-
C ific to ground state CC3 calculations and works with plain integrals.
C DHBIAJK3 depends on prior disk sorting of W, but does not require W
C to be totally symmetric.
C
C
C     IRREPT --- Overall symmetry of T (and Z).
C     LZOFF  --- Target offset (Z is on LZOFF+1,...,LZOFF+4).
C     LWOFF  --- W offset (W is on LZOFF+1,...,LZOFF+6).
C     LT1    --- T list.
C     LT1OFF --- T offset, i.e. T is on parts LT1OFF+1,LT1OFF+2 of LT1.
C
C     Important assumptions :
C
C     Z lists have orders jk;lc, apart from third, which is JkCl.
C     Only fourth Z list is computed for RHF.
C     W lists follow same patterm as 54-59 and have same sign quirk,
C     i.e. 1,2,5,6 must be negated.
C     IRREPW is assumed to be 1 in this routine.
C
C     The skips to 95,195, etc are very important.
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
      DOUBLE PRECISION CORE
      INTEGER MAXCOR,IUHF,IRREPT,LZOFF,LWOFF,LT1,LT1OFF
C-----------------------------------------------------------------------
      INTEGER IINTLN,IFLTLN,IINTFP,IALONE,IBITWD,
     &        POP,VRT,NT,NFEA,NFMI,NSTART,NIRREP,IRREPY,DIRPRD,
     &        IRPDPD,ISYTYP,ID,ISYMOFF
C-----------------------------------------------------------------------
      DOUBLE PRECISION ONE,ONEM,ZILCH
      INTEGER LISTZ,LISTW
      INTEGER I0T1
      INTEGER C,CL,L,LC,ISPIN,DISSIZW,NDISW
      INTEGER I000,I010,I020,I030,NEED,ISIZW,IOFFW,IOFFT,IOFFZ,IOFFW0
      INTEGER IRPJ,IRPK,IRPJK,IRPF,IRPC,IRPL,IRPCL,IRPLC,IRPFK,IRPFJ
      INTEGER ISYMSZ
C-----------------------------------------------------------------------
      DIMENSION CORE(1)
      DIMENSION I0T1(2)
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
C
C-----------------------------------------------------------------------
C     Read T1 amplitudes.
C-----------------------------------------------------------------------
C
      I0T1(1) = 1
      I0T1(2) = I0T1(1) + IRPDPD(IRREPT, 9)
      CALL GETLST(CORE(I0T1(1)),1,1,1,LT1OFF + 1       ,LT1)
      CALL GETLST(CORE(I0T1(2)),1,1,1,LT1OFF + 1 + IUHF,LT1)
C
C-----------------------------------------------------------------------
C     Z(Lc,Jk).
C-----------------------------------------------------------------------
      LISTZ  = LZOFF + 4
C
C     Z(Lc,Jk) =  t(F,J) * W(LcFk). W ordering is FL;ck.
C
      LISTW  =  LWOFF + 3
      ISIZW = ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
C
      I000    = I0T1(2) + IRPDPD(IRREPT,10)
      I010    = I000    + ISIZW
      I020    = I010    + ISIZW
C     Note I020 length can be cut.
      I030    = I020    + ISIZW
      NEED    = I030 * IINTFP
      IF(NEED.GT.MAXCOR) CALL INSMEM('CC3W4RING',NEED,MAXCOR)
C
C --- Pick up W(LcFk), ordered (FL,ck). Reorder to Fk;cL ---
C
      CALL GETALL(CORE(I010),ISIZW,1,LISTW)
      CALL SSTGEN(CORE(I010),CORE(I000),ISIZW,VRT(1,1),POP(1,1),
     &            VRT(1,2),POP(1,2),CORE(I020),1,'1432')
C
      IOFFW0 = I000
      DO 100 IRPCL=1,NIRREP
C
      IRPLC = IRPCL
      IRPFK = IRPCL
      IRPJK = DIRPRD(IRREPT,IRPCL)
C
      DISSIZW = IRPDPD(IRPFK,11)
      NDISW   = IRPDPD(IRPCL,12)
C
      IF(DISSIZW.EQ.0.OR.NDISW.EQ.0) GOTO 100
      IF(IRPDPD(IRPJK,14).EQ.0)      GOTO  95
C
C     I0T1 --- T1(F,J)
C     I000 --- W(LcFk) ordered (Fk,cL)
C     I010 --- Z(J,k)
C     I020 --- Current Hbar element, read from disk, dist. size is J,k
C
      I020 = I010 + IRPDPD(IRPJK,14)
      I030 = I020 + IRPDPD(IRPJK,14)
C
      DO  90  IRPL=1,NIRREP
      IRPC  = DIRPRD(IRPCL,IRPL)
C
      IF(POP(IRPL,1).EQ.0.OR.VRT(IRPC,2).EQ.0) GOTO 90
C
      DO  80     L=1,POP(IRPL,1)
      DO  70     C=1,VRT(IRPC,2)
C
      CL = ISYMOFF(IRPL,IRPCL,12) + (L-1)*VRT(IRPC,2) + C - 1
C
      CALL ZERO(CORE(I010),IRPDPD(IRPJK,14))
C
      DO  60  IRPF=1,NIRREP
      IRPJ  = DIRPRD(IRREPT,IRPF)
      IRPK  = DIRPRD(IRPJ,IRPJK)
C
      IF(VRT(IRPF,1).EQ.0.OR.POP(IRPJ,1).EQ.0.OR.POP(IRPK,2).EQ.0) 
     &                                                     GOTO 60
      IOFFT = ISYMOFF(IRPJ,IRREPT, 9) + I0T1(1) - 1
      IOFFW = IOFFW0 + (CL-1)*DISSIZW +
     &        ISYMOFF(IRPK,IRPFK,11)           - 1
      IOFFZ = ISYMOFF(IRPK,IRPJK,14) + I010    - 1
C
C     T(F,J) * W(F,k)
C
      CALL XGEMM('T','N',POP(IRPJ,1),POP(IRPK,2),
     &                               VRT(IRPF,1),ONE,
     &           CORE(IOFFT),VRT(IRPF,1),
     &           CORE(IOFFW),VRT(IRPF,1),ZILCH,
     &           CORE(IOFFZ),POP(IRPJ,1))
C
   60 CONTINUE
C
      IRPLC = IRPCL
      LC = ISYMOFF(IRPC,IRPLC,18) + (C-1)*POP(IRPL,1) + L - 1
      CALL GETLST(CORE(I020),LC,1,2,IRPLC,LISTZ)
      CALL  SAXPY(IRPDPD(IRPJK,14),ONE,CORE(I010),1,CORE(I020),1)
      CALL PUTLST(CORE(I020),LC,1,2,IRPLC,LISTZ)
C
   70 CONTINUE
   80 CONTINUE
   90 CONTINUE
   95 CONTINUE
      IOFFW0 = IOFFW0 + DISSIZW * NDISW
  100 CONTINUE
C
C     Z(Lc,Jk) =  t(f,k) * <Lc||Jf>
C     UHF => List 26. Ordering is fL;cJ. Reorder to fJ;cL
C     RHF => List 25. As above.
C
C     Z(Lc,Jk) = - t(f,k) * W(LcfJ). Order is fL;cJ. Reorder to fJ;cL.
C
      LISTW  = LWOFF + 5 + IUHF
      ISIZW = ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
C
      I000    = I0T1(2) + IRPDPD(IRREPT,10)
      I010    = I000    + ISIZW
      I020    = I010    + ISIZW
C     Note I020 length can be cut.
      I030    = I020    + ISIZW
      NEED    = I030 * IINTFP
      IF(NEED.GT.MAXCOR) CALL INSMEM('GCC3W4RING',NEED,MAXCOR)
C
      CALL GETALL(CORE(I010),ISIZW,1,LISTW)
      CALL SSTGEN(CORE(I010),CORE(I000),ISIZW,VRT(1,2),POP(1,1),
     &            VRT(1,2),POP(1,1),CORE(I020),1,'1432')
C
      IOFFW0 = I000
      DO 200 IRPCL=1,NIRREP
C
      IRPLC = IRPCL
      IRPFJ = IRPCL
      IRPJK = DIRPRD(IRREPT,IRPCL)
C
      DISSIZW = IRPDPD(IRPFJ,12)
      NDISW   = IRPDPD(IRPCL,12)
C
      IF(DISSIZW.EQ.0.OR.NDISW.EQ.0) GOTO 200
      IF(IRPDPD(IRPJK,14).EQ.0)      GOTO 195
C
C     I0T1 --- T1(f,k)
C     I000 --- W(LcfJ) ordered (fJ,cL)
C     I010 --- Z(J,k)
C     I020 --- Current Hbar element, read from disk, dist. size is J,k
C
      I020 = I010 + IRPDPD(IRPJK,14)
      I030 = I020 + IRPDPD(IRPJK,14)
C
      DO 190  IRPL=1,NIRREP
      IRPC  = DIRPRD(IRPCL,IRPL)
C
      IF(POP(IRPL,1).EQ.0.OR.VRT(IRPC,2).EQ.0) GOTO 190
C
      DO 180     L=1,POP(IRPL,1)
      DO 170     C=1,VRT(IRPC,2)
C
      CL = ISYMOFF(IRPL,IRPCL,12) + (L-1)*VRT(IRPC,2) + C - 1
C
      CALL ZERO(CORE(I010),IRPDPD(IRPJK,14))
C
      DO 160  IRPF=1,NIRREP
      IRPK  = DIRPRD(IRREPT,IRPF)
      IRPJ  = DIRPRD(IRPK,IRPJK)
C
      IF(VRT(IRPF,2).EQ.0.OR.POP(IRPK,2).EQ.0.OR.POP(IRPJ,1).EQ.0) 
     &                                                    GOTO 160
C
      IOFFT = ISYMOFF(IRPK,IRREPT,10) + I0T1(2) - 1
      IOFFW = IOFFW0 + (CL-1)*DISSIZW +
     &        ISYMOFF(IRPJ,IRPFJ,12)           - 1
      IOFFZ = ISYMOFF(IRPK,IRPJK,14) + I010    - 1
C
C     - T(f,k) * W(f,J). Plus since "58,59" lists have minus sign.
C
      CALL XGEMM('T','N',POP(IRPJ,1),POP(IRPK,2),
     &                               VRT(IRPF,2),ONE,
     &           CORE(IOFFW),VRT(IRPF,2),
     &           CORE(IOFFT),VRT(IRPF,2),ZILCH,
     &           CORE(IOFFZ),POP(IRPJ,1))
C
  160 CONTINUE
C
      IRPLC = IRPCL
      LC = ISYMOFF(IRPC,IRPLC,18) + (C-1)*POP(IRPL,1) + L - 1
      CALL GETLST(CORE(I020),LC,1,2,IRPLC,LISTZ)
      CALL  SAXPY(IRPDPD(IRPJK,14),ONE,CORE(I010),1,CORE(I020),1)
      CALL PUTLST(CORE(I020),LC,1,2,IRPLC,LISTZ)
C
  170 CONTINUE
  180 CONTINUE
  190 CONTINUE
  195 CONTINUE
      IOFFW0 = IOFFW0 + DISSIZW * NDISW
  200 CONTINUE
C
      IF(IUHF.EQ.0) RETURN
C
C-----------------------------------------------------------------------
C     Z(Cl,Jk).
C-----------------------------------------------------------------------
C
      LISTZ = LZOFF + 3
C
C     Z(lC,Jk) = t(F,J) * W(lCFk). Ordering is Fl;Ck. Reorder to Fk;Cl.
C     Z(Cl,Jk) = - t(F,J) * W(lCFk).
C
      LISTW  = LWOFF + 5
      ISIZW = ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
C
      I000    = I0T1(2) + IRPDPD(IRREPT,10)
      I010    = I000    + ISIZW
      I020    = I010    + ISIZW
C     Note I020 length can be cut.
      I030    = I020    + ISIZW
      NEED    = I030 * IINTFP
      IF(NEED.GT.MAXCOR) CALL INSMEM('GCC3W4RING',NEED,MAXCOR)
C
C --- Pick up <lC|kF> ordered (Fl;Ck). Reorder to (Fk;Cl). ---
C
      CALL GETALL(CORE(I010),ISIZW,1,LISTW)
      CALL SSTGEN(CORE(I010),CORE(I000),ISIZW,VRT(1,1),POP(1,2),
     &            VRT(1,1),POP(1,2),CORE(I020),1,'1432')
C
      IOFFW0 = I000
      DO 300 IRPCL=1,NIRREP
C
      IRPFK = IRPCL
      IRPJK = DIRPRD(IRREPT,IRPCL)
C
      DISSIZW = IRPDPD(IRPFK,11)
      NDISW   = IRPDPD(IRPCL,11)
C
      IF(DISSIZW.EQ.0.OR.NDISW.EQ.0) GOTO 300
      IF(IRPDPD(IRPJK,14).EQ.0)      GOTO 295
C
C     I0T1 --- T1(F,J)
C     I000 --- W(lCFk) ordered (Fk,Cl)
C     I010 --- Z(J,k)
C     I020 --- Current Hbar element, read from disk, dist. size is J,k
C
      I020 = I010 + IRPDPD(IRPJK,14)
      I030 = I020 + IRPDPD(IRPJK,14)
C
      DO 290  IRPL=1,NIRREP
      IRPC  = DIRPRD(IRPCL,IRPL)
C
      IF(POP(IRPL,2).EQ.0.OR.VRT(IRPC,1).EQ.0) GOTO 290
C
      DO 280     L=1,POP(IRPL,2)
      DO 270     C=1,VRT(IRPC,1)
C
      CL = ISYMOFF(IRPL,IRPCL,11) + (L-1)*VRT(IRPC,1) + C - 1
C
      CALL ZERO(CORE(I010),IRPDPD(IRPJK,14))
C
      DO 260  IRPF=1,NIRREP
      IRPJ  = DIRPRD(IRREPT,IRPF)
      IRPK  = DIRPRD(IRPJ,IRPJK)
C
      IF(VRT(IRPF,1).EQ.0.OR.POP(IRPJ,1).EQ.0.OR.POP(IRPK,2).EQ.0) 
     &                                                     GOTO 260
C
      IOFFT = ISYMOFF(IRPJ,IRREPT, 9) + I0T1(1) - 1
      IOFFW = IOFFW0 + (CL-1)*DISSIZW +
     &        ISYMOFF(IRPK,IRPFK,11)           - 1
      IOFFZ = ISYMOFF(IRPK,IRPJK,14) + I010    - 1
C
C     T(F,J) * W(F,k)
C
      CALL XGEMM('T','N',POP(IRPJ,1),POP(IRPK,2),
     &                               VRT(IRPF,1),ONE,
     &           CORE(IOFFT),VRT(IRPF,1),
     &           CORE(IOFFW),VRT(IRPF,1),ZILCH,
     &           CORE(IOFFZ),POP(IRPJ,1))
C
  260 CONTINUE
C
      CALL GETLST(CORE(I020),CL,1,2,IRPCL,LISTZ)
      CALL  SAXPY(IRPDPD(IRPJK,14),ONE,CORE(I010),1,CORE(I020),1)
      CALL PUTLST(CORE(I020),CL,1,2,IRPCL,LISTZ)
C
  270 CONTINUE
  280 CONTINUE
  290 CONTINUE
  295 CONTINUE
      IOFFW0 = IOFFW0 + DISSIZW * NDISW
  300 CONTINUE
C
C     Z(lC,Jk) =  t(f,k) * <lC||Jf> = -t(f,k) * <Cf|Jl>
C     List 22. Ordering is fJ;Cl.
C     Z(lC,Jk) = - t(f,k) * W(lCfJ). Order is fl;CJ. Reorder to fJ;Cl.
C
      LISTW = LWOFF + 4
      ISIZW = ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
C
      I000    = I0T1(2) + IRPDPD(IRREPT,10)
      I010    = I000    + ISIZW
      I020    = I010    + ISIZW
C     Note I020 length can be cut.
      I030    = I020    + ISIZW
      NEED    = I030 * IINTFP
      IF(NEED.GT.MAXCOR) CALL INSMEM('GCC3W4RING',NEED,MAXCOR)
C
      CALL GETALL(CORE(I010),ISIZW,1,LISTW)
      CALL SSTGEN(CORE(I010),CORE(I000),ISIZW,VRT(1,2),POP(1,2),
     &            VRT(1,1),POP(1,1),CORE(I020),1,'1432')
C
      IOFFW0 = I000
      DO 400 IRPCL=1,NIRREP
C
      IRPFJ = IRPCL
      IRPJK = DIRPRD(IRREPT,IRPCL)
C
      DISSIZW = IRPDPD(IRPFJ,12)
      NDISW   = IRPDPD(IRPCL,11)
C
      IF(DISSIZW.EQ.0.OR.NDISW.EQ.0) GOTO 400
      IF(IRPDPD(IRPJK,14).EQ.0)      GOTO 395
C
C     I0T1 --- T1(f,k)
C     I000 --- W(lCJf) ordered (fJ,Cl)
C     I010 --- Z(J,k)
C     I020 --- Current Hbar element, read from disk, dist. size is J,k
C
      I020 = I010 + IRPDPD(IRPJK,14)
      I030 = I020 + IRPDPD(IRPJK,14)
C
      DO 390  IRPL=1,NIRREP
      IRPC  = DIRPRD(IRPCL,IRPL)
C
      IF(POP(IRPL,2).EQ.0.OR.VRT(IRPC,1).EQ.0) GOTO 390
C
      DO 380     L=1,POP(IRPL,2)
      DO 370     C=1,VRT(IRPC,1)
C
      CL = ISYMOFF(IRPL,IRPCL,11) + (L-1)*VRT(IRPC,1) + C - 1
C
      CALL ZERO(CORE(I010),IRPDPD(IRPJK,14))
C
      DO 360  IRPF=1,NIRREP
      IRPK  = DIRPRD(IRREPT,IRPF)
      IRPJ  = DIRPRD(IRPK,IRPJK)
C
      IF(VRT(IRPF,2).EQ.0.OR.POP(IRPJ,1).EQ.0.OR.POP(IRPK,2).EQ.0) 
     &                                                     GOTO 360
C
      IOFFT = ISYMOFF(IRPK,IRREPT,10) + I0T1(2) - 1
      IOFFW = IOFFW0 + (CL-1)*DISSIZW +
     &        ISYMOFF(IRPJ,IRPFJ,12)           - 1
      IOFFZ = ISYMOFF(IRPK,IRPJK,14) + I010    - 1
C
C     T(f,k) * W(f,J)
C
      CALL XGEMM('T','N',POP(IRPJ,1),POP(IRPK,2),
     &                               VRT(IRPF,2),ONE,
     &           CORE(IOFFW),VRT(IRPF,2),
     &           CORE(IOFFT),VRT(IRPF,2),ZILCH,
     &           CORE(IOFFZ),POP(IRPJ,1))
C
  360 CONTINUE
C
      CALL GETLST(CORE(I020),CL,1,2,IRPCL,LISTZ)
      CALL  SAXPY(IRPDPD(IRPJK,14),ONE,CORE(I010),1,CORE(I020),1)
      CALL PUTLST(CORE(I020),CL,1,2,IRPCL,LISTZ)
C
  370 CONTINUE
  380 CONTINUE
  390 CONTINUE
  395 CONTINUE
      IOFFW0 = IOFFW0 + DISSIZW * NDISW
  400 CONTINUE
C
C
C-----------------------------------------------------------------------
C     Z(LC,JK), Z(lc,jk).
C-----------------------------------------------------------------------
      DO 600 ISPIN=1,2
C
      LISTZ = LZOFF + ISPIN
C
C     Z(LC,JK) = -t(F,J) * <LC||KF>
C     => List 22 + ISPIN. Ordering is FL;CK. Reorder to FK;CL.
C     Permutation and compression handled by ASSYM2.
C     Z(LC,JK) = t(F,J) * W(LCFK). Order is FL;CK. Reorder to FK;CL.
C
      LISTW  = LWOFF + ISPIN
      ISIZW = ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
C
      I000    = I0T1(2) + IRPDPD(IRREPT,10)
      I010    = I000    + ISIZW
      I020    = I010    + ISIZW
C     Note I020 length can be cut.
      I030    = I020    + ISIZW
      NEED    = I030 * IINTFP
      IF(NEED.GT.MAXCOR) CALL INSMEM('GCC3W4RING',NEED,MAXCOR)
C
C --- Pick up <LC|KF> ordered (FL;CK). Reorder to (FK;CL). ---
C
      CALL GETALL(CORE(I010),ISIZW,1,LISTW)
      CALL SSTGEN(CORE(I010),CORE(I000),ISIZW,
     &            VRT(1,ISPIN),POP(1,ISPIN),
     &            VRT(1,ISPIN),POP(1,ISPIN),CORE(I020),1,'1432')
C
      IOFFW0 = I000
      DO 500 IRPCL=1,NIRREP
C
      IRPLC = IRPCL
      IRPFK = IRPCL
      IRPJK = DIRPRD(IRREPT,IRPCL)
C
      DISSIZW = IRPDPD(IRPFK,8+ISPIN)
      NDISW   = IRPDPD(IRPCL,8+ISPIN)
C
      IF(DISSIZW.EQ.0.OR.NDISW.EQ.0)  GOTO 500
      IF(IRPDPD(IRPJK, 2+ISPIN).EQ.0) GOTO 495
      IF(IRPDPD(IRPJK,20+ISPIN).EQ.0) GOTO 495
C
C     I0T1 --- T1(F,J)
C     I000 --- W(LCKF) ordered (FK,CL)
C     I010 --- Z(J,K)
C     I020 --- Current Hbar element, read from disk, dist. size is J<K
C
      I020 = I010 + IRPDPD(IRPJK,20+ISPIN)
      I030 = I020 + IRPDPD(IRPJK,20+ISPIN)
C
      DO 490  IRPL=1,NIRREP
      IRPC  = DIRPRD(IRPCL,IRPL)
C
      IF(POP(IRPL,ISPIN).EQ.0.OR.VRT(IRPC,ISPIN).EQ.0) GOTO 490
C
      DO 480     L=1,POP(IRPL,ISPIN)
      DO 470     C=1,VRT(IRPC,ISPIN)
C
      CL = ISYMOFF(IRPL,IRPCL,8+ISPIN) + (L-1)*VRT(IRPC,ISPIN) + C - 1
C
      DO 460  IRPF=1,NIRREP
      IRPJ  = DIRPRD(IRREPT,IRPF)
      IRPK  = DIRPRD(IRPJ,IRPJK)
C
      IOFFT = ISYMOFF(IRPJ,IRREPT, 8+ISPIN) + I0T1(ISPIN) - 1
      IOFFW = IOFFW0 + (CL-1)*DISSIZW +
     &        ISYMOFF(IRPK,IRPFK, 8+ISPIN)           - 1
      IOFFZ = ISYMOFF(IRPK,IRPJK,20+ISPIN) + I010    - 1
C
C     T(F,J) * W(F,K)
C
      CALL XGEMM('T','N',POP(IRPJ,ISPIN),POP(IRPK,ISPIN),
     &                                   VRT(IRPF,ISPIN),ONEM,
     &           CORE(IOFFT),VRT(IRPF,ISPIN),
     &           CORE(IOFFW),VRT(IRPF,ISPIN),ZILCH,
     &           CORE(IOFFZ),POP(IRPJ,ISPIN))
C
  460 CONTINUE
C
      CALL ASSYM2(IRPJK,POP(1,ISPIN),1,CORE(I010))
C
      IRPLC = IRPCL
      LC = ISYMOFF(IRPC,IRPLC,15+ISPIN) + 
     &     (C-1)*POP(IRPL,ISPIN) + L - 1
      CALL GETLST(CORE(I020),LC,1,2,IRPLC,LISTZ)
      CALL  SAXPY(IRPDPD(IRPJK, 2+ISPIN),ONE,CORE(I010),1,
     &                                       CORE(I020),1)
      CALL PUTLST(CORE(I020),LC,1,2,IRPLC,LISTZ)
C
  470 CONTINUE
  480 CONTINUE
  490 CONTINUE
  495 CONTINUE
      IOFFW0 = IOFFW0 + DISSIZW * NDISW
  500 CONTINUE
  600 CONTINUE
      RETURN
      END
