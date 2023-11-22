      SUBROUTINE CC3W4RING(CORE,MAXCOR,IUHF)
C
C This subroutine is used to compute the contribution
C
C     Z(lc,jk) = - \sum_f [ t(f,j) * <lc||kf> + t(f,k) * <lc||jf> ]
C     Z        =            T * W
C
C to the Hbar element. It is called in CC3 calculations. A similar
C contraction is performed by W4RING, but CC3W4RING needs only integrals.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POP,VRT,DIRPRD,C,CL,DISSIZW
      DIMENSION CORE(1),I0T1(2),LENT1(8),IOFFT1(8,2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFEA(2),NFMI(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
C
      DATA ONE,ONEM,ZILCH /1.0D+00,-1.0D+00,0.0D+00/
C
C --- Read T1 amplitudes ---
C
      LISTT1 =  90
      I0T1(1) = 1
      I0T1(2) = I0T1(1) + NT(1)
      CALL GETLST(CORE(I0T1(1)),1,1,1,1     ,LISTT1)
      CALL GETLST(CORE(I0T1(2)),1,1,1,1+IUHF,LISTT1)
C
C --- Make offset arrays for T1 vectors ---
C
      DO 20 ISPIN=1,2
      DO 10 IRPJ =1,NIRREP
        LENT1(IRPJ) = VRT(IRPJ,ISPIN) * POP(IRPJ,ISPIN)
        IF(IRPJ.EQ.1)THEN
          IOFFT1(IRPJ,ISPIN) = 0
        ELSE
          IOFFT1(IRPJ,ISPIN) = IOFFT1(IRPJ-1,ISPIN) + LENT1(IRPJ-1)
        ENDIF
   10 CONTINUE
   20 CONTINUE
C
      LISTZ  = 110
C
C     Z(Lc,Jk) = -t(F,J) * <Lc||kF> = t(F,J) * (LF|kc)
C              =  t(F,J) * <Fc| Lk>
C     => List 21. Ordering is Fk;cL
C
      LISTW  =  21
      ISIZW = ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
C
      I000    = I0T1(2) + NT(2)
      I010    = I000    + ISIZW
      I020    = I010    + ISIZW
C     Note I020 length can be cut.
      I030    = I020    + ISIZW
      NEED    = I030 * IINTFP
      IF(NEED.GT.MAXCOR) CALL INSMEM('CC3W4RING',NEED,MAXCOR)
C
C --- Pick up <Fc|Lk> ordered (Fk,cL) ---
C
      CALL GETALL(CORE(I000),ISIZW,1,LISTW)
C
      IOFFW0 = I000
      DO 100 IRPCL=1,NIRREP
C
      IRPLC = IRPCL
      IRPFK = IRPCL
      IRPJK = IRPCL
C
      DISSIZW = IRPDPD(IRPFK,11)
      NDISW   = IRPDPD(IRPCL,12)
C
      IF(DISSIZW.EQ.0.OR.NDISW.EQ.0) GOTO 100
      IF(IRPDPD(IRPJK,14).EQ.0)      GOTO 100
C
C     I0T1 --- T1(F,J)
C     I000 --- <Lc||kF> ordered (Fk,cL)
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
      DO  60  IRPF=1,NIRREP
      IRPJ  = IRPF
      IRPK  = DIRPRD(IRPJ,IRPJK)
C
      IOFFT = ISYMOFF(IRPJ,    1, 9) + I0T1(1) - 1
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
      IOFFW0 = IOFFW0 + DISSIZW * NDISW
  100 CONTINUE
C
C     Z(Lc,Jk) =  t(f,k) * <Lc||Jf>
C     UHF => List 26. Ordering is fL;cJ. Reorder to fJ;cL
C     RHF => List 25. As above.
C
      LISTW  =  25 + IUHF
      ISIZW = ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
C
      I000    = I0T1(2) + NT(2)
      I010    = I000    + ISIZW
      I020    = I010    + ISIZW
C     Note I020 length can be cut.
      I030    = I020    + ISIZW
      NEED    = I030 * IINTFP
      IF(NEED.GT.MAXCOR) CALL INSMEM('CC3W4RING',NEED,MAXCOR)
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
      IRPJK = IRPCL
C
      DISSIZW = IRPDPD(IRPFJ,12)
      NDISW   = IRPDPD(IRPCL,12)
C
      IF(DISSIZW.EQ.0.OR.NDISW.EQ.0) GOTO 200
      IF(IRPDPD(IRPJK,14).EQ.0)      GOTO 200
C
C     I0T1 --- T1(f,k)
C     I000 --- <Lc||Jf> ordered (FJ,cL)
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
      DO 160  IRPF=1,NIRREP
      IRPK  = IRPF
      IRPJ  = DIRPRD(IRPK,IRPJK)
C
      IOFFT = ISYMOFF(IRPK,    1,10) + I0T1(2) - 1
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
      IOFFW0 = IOFFW0 + DISSIZW * NDISW
  200 CONTINUE
C
      IF(IUHF.EQ.0) RETURN
C
      LISTZ = 109
C
C     Z(lC,Jk) = -t(F,J) * <lC||kF>
C     => List 25. Ordering is Fl;Ck. Reorder to Fk;Cl.
C
      LISTW  =  25
      ISIZW = ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
C
      I000    = I0T1(2) + NT(2)
      I010    = I000    + ISIZW
      I020    = I010    + ISIZW
C     Note I020 length can be cut.
      I030    = I020    + ISIZW
      NEED    = I030 * IINTFP
      IF(NEED.GT.MAXCOR) CALL INSMEM('CC3W4RING',NEED,MAXCOR)
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
      IRPLC = IRPCL
      IRPFK = IRPCL
      IRPJK = IRPCL
C
      DISSIZW = IRPDPD(IRPFK,11)
      NDISW   = IRPDPD(IRPCL,11)
C
      IF(DISSIZW.EQ.0.OR.NDISW.EQ.0) GOTO 300
      IF(IRPDPD(IRPJK,14).EQ.0)      GOTO 300
C
C     I0T1 --- T1(F,J)
C     I000 --- <lC||kF> ordered (Fk,Cl)
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
C     CL = ISYMOFF(IRPL,IRPCL,12) + (L-1)*VRT(IRPC,2) + C - 1
      CL = ISYMOFF(IRPL,IRPCL,11) + (L-1)*VRT(IRPC,1) + C - 1
C
      DO 260  IRPF=1,NIRREP
      IRPJ  = IRPF
      IRPK  = DIRPRD(IRPJ,IRPJK)
C
      IOFFT = ISYMOFF(IRPJ,    1, 9) + I0T1(1) - 1
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
C     IRPLC = IRPCL
C     LC = ISYMOFF(IRPC,IRPLC,18) + (C-1)*POP(IRPL,1) + L - 1
      CALL GETLST(CORE(I020),CL,1,2,IRPLC,LISTZ)
      CALL  SAXPY(IRPDPD(IRPJK,14),ONE,CORE(I010),1,CORE(I020),1)
      CALL PUTLST(CORE(I020),CL,1,2,IRPLC,LISTZ)
C
  270 CONTINUE
  280 CONTINUE
  290 CONTINUE
      IOFFW0 = IOFFW0 + DISSIZW * NDISW
  300 CONTINUE
C
C     Z(lC,Jk) =  t(f,k) * <lC||Jf> = -t(f,k) * <Cf|Jl>
C     List 22. Ordering is fJ;Cl.
C
      LISTW  =  22
      ISIZW = ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
C
      I000    = I0T1(2) + NT(2)
      I010    = I000    + ISIZW
      I020    = I010    + ISIZW
C     Note I020 length can be cut.
      I030    = I020    + ISIZW
      NEED    = I030 * IINTFP
      IF(NEED.GT.MAXCOR) CALL INSMEM('CC3W4RING',NEED,MAXCOR)
C
      CALL GETALL(CORE(I000),ISIZW,1,LISTW)
C
      IOFFW0 = I000
      DO 400 IRPCL=1,NIRREP
C
      IRPLC = IRPCL
      IRPFJ = IRPCL
      IRPJK = IRPCL
C
      DISSIZW = IRPDPD(IRPFJ,12)
      NDISW   = IRPDPD(IRPCL,11)
C
      IF(DISSIZW.EQ.0.OR.NDISW.EQ.0) GOTO 400
      IF(IRPDPD(IRPJK,14).EQ.0)      GOTO 400
C
C     I0T1 --- T1(f,k)
C     I000 --- <lC||Jf> ordered (fJ,Cl)
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
      DO 360  IRPF=1,NIRREP
      IRPK  = IRPF
      IRPJ  = DIRPRD(IRPK,IRPJK)
C
      IOFFT = ISYMOFF(IRPK,    1,10) + I0T1(2) - 1
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
C     IRPLC = IRPCL
C     LC = ISYMOFF(IRPC,IRPLC,18) + (C-1)*POP(IRPL,1) + L - 1
      CALL GETLST(CORE(I020),CL,1,2,IRPLC,LISTZ)
      CALL  SAXPY(IRPDPD(IRPJK,14),ONE,CORE(I010),1,CORE(I020),1)
      CALL PUTLST(CORE(I020),CL,1,2,IRPLC,LISTZ)
C
  370 CONTINUE
  380 CONTINUE
  390 CONTINUE
      IOFFW0 = IOFFW0 + DISSIZW * NDISW
  400 CONTINUE
C
      DO 600 ISPIN=1,2
C
      LISTZ = 106 + ISPIN
C
C     Z(LC,JK) = -t(F,J) * <LC||KF>
C     => List 22 + ISPIN. Ordering is FL;CK. Reorder to FK;CL.
C     Permutation and compression handled by ASSYM2.
C
      LISTW  =  22 + ISPIN
      ISIZW = ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
C
      I000    = I0T1(2) + NT(2)
      I010    = I000    + ISIZW
      I020    = I010    + ISIZW
C     Note I020 length can be cut.
      I030    = I020    + ISIZW
      NEED    = I030 * IINTFP
      IF(NEED.GT.MAXCOR) CALL INSMEM('CC3W4RING',NEED,MAXCOR)
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
      IRPJK = IRPCL
C
      DISSIZW = IRPDPD(IRPFK,8+ISPIN)
      NDISW   = IRPDPD(IRPCL,8+ISPIN)
C
      IF(DISSIZW.EQ.0.OR.NDISW.EQ.0)  GOTO 500
      IF(IRPDPD(IRPJK, 2+ISPIN).EQ.0) GOTO 500
      IF(IRPDPD(IRPJK,20+ISPIN).EQ.0) GOTO 500
C
C     I0T1 --- T1(F,J)
C     I000 --- <LC||KF> ordered (FK,CL)
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
      IRPJ  = IRPF
      IRPK  = DIRPRD(IRPJ,IRPJK)
C
      IOFFT = ISYMOFF(IRPJ,    1, 8+ISPIN) + I0T1(ISPIN) - 1
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
      IOFFW0 = IOFFW0 + DISSIZW * NDISW
  500 CONTINUE
  600 CONTINUE
C
      RETURN
      END
