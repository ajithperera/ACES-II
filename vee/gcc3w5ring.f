      SUBROUTINE GCC3W5RING(CORE,MAXCOR,IUHF,IRREPT,
     &                      LZOFF,LWOFF,LT1,LT1OFF)
C
C This subroutine is used to compute the contribution
C
C     Z(bc,dk) =   \sum_o [ t(c,o) * W(obdk) + t(b,o) * W(ocdk) ]
C     Z        =            T * W
C
C IRREPW is assumed to be 1. W(obdk) is assumed to be stored d,o; b,k.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POP,VRT,DIRPRD,D,DK
      DIMENSION CORE(1),I0T1(2)
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
      I0T1(1) = 1
      I0T1(2) = I0T1(1) + IRPDPD(IRREPT, 9)
      CALL GETLST(CORE(I0T1(1)),1,1,1,LT1OFF + 1     ,LT1)
      CALL GETLST(CORE(I0T1(2)),1,1,1,LT1OFF + 1+IUHF,LT1)
C
C     Z(Bc,Dk) = t(c,n) * W(nBDk) - t(B,N) * W(NcDk)
C
      LISTZ  = LZOFF + 4
C
C     I think the following line is wrong for UHF. We should ALWAYS pick
C     up list LWOFF + 5.
C     LISTW  = LWOFF + 5 + IUHF
      LISTW  = LWOFF + 5
      ISIZW = ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
C
      I000    = I0T1(2) + IRPDPD(IRREPT,10)
      I010    = I000    + ISIZW
      I020    = I010    + ISIZW
C     Note I020 length can be cut.
      I030    = I020    + ISIZW
      NEED    = I030 * IINTFP
      IF(NEED.GT.MAXCOR) CALL INSMEM('GCC3W5RING',NEED,MAXCOR)
C
C --- Pick up W(nBDk) ordered (Dn,Bk). Reorder to (Bn,Dk) ---
C
      CALL GETALL(CORE(I010),ISIZW,1,LISTW)
      CALL SSTGEN(CORE(I010),CORE(I000),ISIZW,VRT(1,1),POP(1,2),
     &            VRT(1,1),POP(1,2),CORE(I020),1,'3214')
C
      IOFFW0 = I000
C
      DO  50 IRPDK=1,NIRREP
C
      IRPBC = DIRPRD(IRREPT,IRPDK)
C
      IF(IRPDPD(IRPDK,11).EQ.0) GOTO 50
C
      IOFFW = IOFFW0
C
      DO  40  IRPK=1,NIRREP
      IRPD  = DIRPRD(IRPDK,IRPK)
C
      IF(POP(IRPK,2).EQ.0.OR.VRT(IRPD,1).EQ.0) GOTO 40
C
      DO  30     K=1,POP(IRPK,2)
      DO  20     D=1,VRT(IRPD,1)
C
      IOFFT = I0T1(2)
C
      I020 = I010 + IRPDPD(IRPBC,13)
      I030 = I020 + IRPDPD(IRPBC,13)
C
C     I0T1 --- T1(c,n)
C     I000 --- W(nBDk) ordered (bn,dk)
C     I010 --- Z(B,c)
C     I020 --- Current Hbar element, read from disk, dist. size is B,c
C
      CALL ZERO(CORE(I010),IRPDPD(IRPBC,13))
      CALL ZERO(CORE(I020),IRPDPD(IRPBC,13))
      CALL ZERO(CORE(I030),IRPDPD(IRPBC,13))
C
      DO  10  IRPN=1,NIRREP
      IRPC  = DIRPRD(IRREPT,IRPN)
      IRPB  = DIRPRD(IRPC,IRPBC)
      IOFFZ = I010 - 1 + ISYMOFF(IRPC,IRPBC,13)
C
C     W(B,n) * T(c,n) (negate since lists 58,59 are minus lists).
C
      IF(VRT(IRPB,1).GT.0.AND.VRT(IRPC,2).GT.0.AND.POP(IRPN,2).GT.0)THEN
      CALL XGEMM('N','T',VRT(IRPB,1),VRT(IRPC,2),
     &                               POP(IRPN,2),ONEM,
     &           CORE(IOFFW),VRT(IRPB,1),
     &           CORE(IOFFT),VRT(IRPC,2),ZILCH,
     &           CORE(IOFFZ),VRT(IRPB,1))
      ENDIF
C
      IOFFT = IOFFT + VRT(IRPC,2) * POP(IRPN,2)
      IOFFW = IOFFW + VRT(IRPB,1) * POP(IRPN,2)
   10 CONTINUE
C
      DK = ISYMOFF(IRPK,IRPDK,11) + (K-1)*VRT(IRPD,1) + D - 1
      IF(IRPDPD(IRPBC,13).GT.0)THEN
      CALL GETLST(CORE(I020),DK,1,2,IRPDK,LISTZ)
      CALL  SAXPY(IRPDPD(IRPBC,13),ONE,CORE(I010),1,CORE(I020),1)
      CALL PUTLST(CORE(I020),DK,1,2,IRPDK,LISTZ)
      ENDIF
C
   20 CONTINUE
   30 CONTINUE
   40 CONTINUE
      IOFFW0 = IOFFW0 + IRPDPD(IRPDK,11) * IRPDPD(IRPDK,11)
   50 CONTINUE
C
C     Z(Bc,Dk) = - t(B,N) * W(NcDk)
C
      LISTZ  = LZOFF + 4
C
      LISTW  = LWOFF + 3
      ISIZW = ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
C
      I000    = I0T1(2) + IRPDPD(IRREPT,10)
      I010    = I000    + ISIZW
      I020    = I010    + ISIZW
C     Note I020 length can be cut.
      I030    = I020    + ISIZW
      NEED    = I030 * IINTFP
      IF(NEED.GT.MAXCOR) CALL INSMEM('GCC3W5RING',NEED,MAXCOR)
C
C --- Pick up W(NcDk) ordered DN;ck. Reorder to cN;Dk.
C
      CALL GETALL(CORE(I010),ISIZW,1,LISTW)
      CALL SSTGEN(CORE(I010),CORE(I000),ISIZW,VRT(1,1),POP(1,1),
     &            VRT(1,2),POP(1,2),CORE(I020),1,'3214')
C
      IOFFW0 = I000
C
      DO 100 IRPDK=1,NIRREP
C
      IRPBC = DIRPRD(IRREPT,IRPDK)
C
      IF(IRPDPD(IRPDK,11).EQ.0) GOTO 100
C
c      DK = 1
C
      IOFFW = IOFFW0
C
      DO  90  IRPK=1,NIRREP
      IRPD  = DIRPRD(IRPDK,IRPK)
C
      IF(POP(IRPK,2).EQ.0.OR.VRT(IRPD,1).EQ.0) GOTO 90
C
      DO  80     K=1,POP(IRPK,2)
      DO  70     D=1,VRT(IRPD,1)
C
      IOFFT = I0T1(1)
C
      I020 = I010 + IRPDPD(IRPBC,13)
      I030 = I020 + IRPDPD(IRPBC,13)
C
C     I0T1 --- T1(B,N)
C     I000 --- W(NcDk) ordered (cN,Dk)
C     I010 --- Z(B,c)
C     I020 --- Current Hbar element, read from disk, dist. size is B,c
C
      CALL ZERO(CORE(I010),IRPDPD(IRPBC,13))
      CALL ZERO(CORE(I020),IRPDPD(IRPBC,13))
      CALL ZERO(CORE(I030),IRPDPD(IRPBC,13))
C
      DO  60  IRPN=1,NIRREP
      IRPB  = DIRPRD(IRREPT,IRPN)
      IRPC  = DIRPRD(IRPB,IRPBC)
      IOFFZ = I010 - 1 + ISYMOFF(IRPB,IRPBC,23)
C
C     - W(c,N) * T(B,N)
C
      IF(VRT(IRPC,2).GT.0.AND.VRT(IRPB,1).GT.0.AND.POP(IRPN,1).GT.0)THEN
      CALL XGEMM('N','T',VRT(IRPC,2),VRT(IRPB,1),
     &                               POP(IRPN,1),ONEM,
     &           CORE(IOFFW),VRT(IRPC,2),
     &           CORE(IOFFT),VRT(IRPB,1),ZILCH,
     &           CORE(IOFFZ),VRT(IRPC,2))
      ENDIF
C
      IOFFT = IOFFT + VRT(IRPB,1) * POP(IRPN,1)
      IOFFW = IOFFW + VRT(IRPC,2) * POP(IRPN,1)
   60 CONTINUE
C
      DK = ISYMOFF(IRPK,IRPDK,11) + (K-1)*VRT(IRPD,1) + D - 1
      IF(IRPDPD(IRPBC,13).GT.0)THEN
      CALL SYMTRA2(IRPBC,VRT(1,2),VRT(1,1),IRPDPD(IRPBC,13),1,
     &             CORE(I010),CORE(I030))
      CALL GETLST(CORE(I020),DK,1,2,IRPDK,LISTZ)
      CALL  SAXPY(IRPDPD(IRPBC,13),ONE,CORE(I030),1,CORE(I020),1)
      CALL PUTLST(CORE(I020),DK,1,2,IRPDK,LISTZ)
      ENDIF
C
c      DK = DK + 1
   70 CONTINUE
   80 CONTINUE
   90 CONTINUE
      IOFFW0 = IOFFW0 + IRPDPD(IRPDK,11)*IRPDPD(IRPDK,11)
  100 CONTINUE
C
      IF(IUHF.EQ.0) RETURN
C
      DO 200 ISPIN=1,2
C
      LISTZ = LZOFF + ISPIN
      LISTW = LWOFF + ISPIN
C
      ISIZW = ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
C
      I000    = I0T1(2) + IRPDPD(IRREPT,10)
      I010    = I000    + ISIZW
      I020    = I010    + ISIZW
C     Note I020 length can be cut.
      I030    = I020    + ISIZW
      NEED    = I030 * IINTFP
      IF(NEED.GT.MAXCOR) CALL INSMEM('CC3W5RING',NEED,MAXCOR)
C
C --- Pick up W(nbdk) ordered (dn,bk). Reorder to (bn,dk) ---
C
      CALL GETALL(CORE(I010),ISIZW,1,LISTW)
      CALL SSTGEN(CORE(I010),CORE(I000),ISIZW,VRT(1,ISPIN),POP(1,ISPIN),
     &            VRT(1,ISPIN),POP(1,ISPIN),CORE(I020),1,'3214')
C
      IOFFW = I000
C
      DO 150 IRPDK=1,NIRREP
C
      IRPBC = DIRPRD(IRREPT,IRPDK)
C
      IF(IRPDPD(IRPDK,8+ISPIN).EQ.0 .OR. IRPDPD(IRPBC,ISPIN).EQ.0)
     &                                                    GOTO 150
C
      DO 140  IRPK=1,NIRREP
      IRPD  = DIRPRD(IRPDK,IRPK)

      DO 130     K=1,POP(IRPK,ISPIN)
      DO 120     D=1,VRT(IRPD,ISPIN)
C
      IOFFT = I0T1(ISPIN)
C
      I020 = I010 + IRPDPD(IRPBC,18+ISPIN)
      I030 = I020 + IRPDPD(IRPBC,   ISPIN)
C
C     I0T1 --- T1(c,n)
C     I000 --- W(nbdk) ordered (bn,dk)
C     I010 --- Z(b,c); after ASSYM2 becomes Z(b<c)
C     I020 --- Current Hbar element, read from disk, dist. size is b<c
C
      CALL ZERO(CORE(I010),IRPDPD(IRPBC,18+ISPIN))
      CALL ZERO(CORE(I020),IRPDPD(IRPBC,   ISPIN))
C
      DO 110  IRPN=1,NIRREP
      IRPC  = DIRPRD(IRREPT,IRPN)
      IRPB  = DIRPRD(IRPC,IRPBC)
C
      IOFFZ = I010 - 1 + ISYMOFF(IRPC,IRPBC,18+ISPIN)
C
C     W(b,n) * T(c,n) (negate because W lists 54,54 are minus).
C
      CALL XGEMM('N','T',VRT(IRPB,ISPIN),VRT(IRPC,ISPIN),
     &                                   POP(IRPN,ISPIN),ONEM,
     &           CORE(IOFFW),VRT(IRPB,ISPIN),
     &           CORE(IOFFT),VRT(IRPC,ISPIN),ZILCH,
     &           CORE(IOFFZ),VRT(IRPB,ISPIN))
C
      IOFFT = IOFFT + VRT(IRPC,ISPIN) * POP(IRPN,ISPIN)
      IOFFW = IOFFW + VRT(IRPB,ISPIN) * POP(IRPN,ISPIN)
  110 CONTINUE
C
      CALL ASSYM2(IRPBC,VRT(1,ISPIN),1,CORE(I010))
      DK = ISYMOFF(IRPK,IRPDK, 8+ISPIN) + (K-1)*VRT(IRPD,ISPIN) + D - 1
      CALL GETLST(CORE(I020),DK,1,2,IRPDK,LISTZ)
      CALL  SAXPY(IRPDPD(IRPBC,ISPIN),ONE,CORE(I010),1,CORE(I020),1)
      CALL PUTLST(CORE(I020),DK,1,2,IRPDK,LISTZ)
C
  120 CONTINUE
  130 CONTINUE
  140 CONTINUE
  150 CONTINUE
  200 CONTINUE
C
C     Z(Bc,dK) =   t(c,n) * W(nBdK) - t(B,N) * W(NcdK)
C     Z(Bc,Kd) = - t(c,n) * W(nBdK) + t(B,N) * W(NcdK)
C
      LISTZ  = LZOFF + 3
C
      LISTW  = LWOFF + 4
      ISIZW = ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
C
      I000    = I0T1(2) + IRPDPD(IRREPT,10)
      I010    = I000    + ISIZW
      I020    = I010    + ISIZW
C     Note I020 length can be cut.
      I030    = I020    + ISIZW
      NEED    = I030 * IINTFP
      IF(NEED.GT.MAXCOR) CALL INSMEM('CC3W5RING',NEED,MAXCOR)
C
C --- Pick up W(nBdK) ordered (dn,BK). Reorder to (Bn,dK) ---
C
      CALL GETALL(CORE(I010),ISIZW,1,LISTW)
      CALL SSTGEN(CORE(I010),CORE(I000),ISIZW,VRT(1,2),POP(1,2),
     &            VRT(1,1),POP(1,1),CORE(I020),1,'3214')
C
      IOFFW = I000
C
      DO 250 IRPDK=1,NIRREP
C
      IRPBC = DIRPRD(IRREPT,IRPDK)
C
      IF(IRPDPD(IRPDK,18).EQ.0 .OR. IRPDPD(IRPBC,13).EQ.0) GOTO 250
C
      DO 240  IRPK=1,NIRREP
      IRPD  = DIRPRD(IRPDK,IRPK)
C
      IF(POP(IRPK,1).EQ.0 .OR. VRT(IRPD,2).EQ.0) GOTO 240
C
      DO 230     K=1,POP(IRPK,1)
      DO 220     D=1,VRT(IRPD,2)
C
      IOFFT = I0T1(2)
C
      I020 = I010 + IRPDPD(IRPBC,13)
      I030 = I020 + IRPDPD(IRPBC,13)
C
C     I0T1 --- T1(c,n)
C     I000 --- W(nBdK) ordered (bn,dk)
C     I010 --- Z(B,c)
C     I020 --- Current Hbar element, read from disk, dist. size is B,c
C
      CALL ZERO(CORE(I010),IRPDPD(IRPBC,13))
      CALL ZERO(CORE(I020),IRPDPD(IRPBC,13))
      CALL ZERO(CORE(I030),IRPDPD(IRPBC,13))
C
      DO 210  IRPN=1,NIRREP
      IRPC  = DIRPRD(IRREPT,IRPN)
      IRPB  = DIRPRD(IRPC,IRPBC)
C
C     - W(B,n) * T(c,n)
C
      IOFFZ = I010 - 1 + ISYMOFF(IRPC,IRPBC,13)
C
      CALL XGEMM('N','T',VRT(IRPB,1),VRT(IRPC,2),
     &                               POP(IRPN,2),ONEM,
     &           CORE(IOFFW),VRT(IRPB,1),
     &           CORE(IOFFT),VRT(IRPC,2),ZILCH,
     &           CORE(IOFFZ),VRT(IRPB,1))
C
      IOFFT = IOFFT + VRT(IRPC,2) * POP(IRPN,2)
      IOFFW = IOFFW + VRT(IRPB,1) * POP(IRPN,2)
c     IOFFZ = IOFFZ + VRT(IRPB,1) * VRT(IRPC,2)
  210 CONTINUE
C
      KD = ISYMOFF(IRPD,IRPDK,18) + K - 1 + (D-1)*POP(IRPK,1)
      CALL GETLST(CORE(I020),KD,1,2,IRPDK,LISTZ)
      CALL  SAXPY(IRPDPD(IRPBC,13),ONE,CORE(I010),1,CORE(I020),1)
      CALL PUTLST(CORE(I020),KD,1,2,IRPDK,LISTZ)
C
  220 CONTINUE
  230 CONTINUE
  240 CONTINUE
  250 CONTINUE
C
C     Z(Bc,Kd) = - t(B,N) * W(NcdK)
C
      LISTZ  = LZOFF + 3
C
      LISTW  = LWOFF + 6
      ISIZW = ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
C
      I000    = I0T1(2) + IRPDPD(IRREPT,10)
      I010    = I000    + ISIZW
      I020    = I010    + ISIZW
C     Note I020 length can be cut.
      I030    = I020    + ISIZW
      NEED    = I030 * IINTFP
      IF(NEED.GT.MAXCOR) CALL INSMEM('CC3W5RING',NEED,MAXCOR)
C
C --- Pick up W(NcdK), ordered (dN,cK). Reorder to (cN,dK).
C
      CALL GETALL(CORE(I010),ISIZW,1,LISTW)
      CALL SSTGEN(CORE(I010),CORE(I000),ISIZW,VRT(1,2),POP(1,1),
     &            VRT(1,2),POP(1,1),CORE(I030),1,'3214')
C
      IOFFW = I000
C
      DO 300 IRPDK=1,NIRREP
C
      IRPBC = DIRPRD(IRREPT,IRPDK)
C
      IF(IRPDPD(IRPDK,18).EQ.0 .OR. IRPDPD(IRPBC,13).EQ.0) GOTO 300
C
      DO 290  IRPK=1,NIRREP
      IRPD  = DIRPRD(IRPDK,IRPK)
C
      IF(POP(IRPK,1).EQ.0 .OR. VRT(IRPD,2).EQ.0) GOTO 290
C
      DO 280     K=1,POP(IRPK,1)
      DO 270     D=1,VRT(IRPD,2)
C
      IOFFT = I0T1(1)
C
      I020 = I010 + IRPDPD(IRPBC,13)
      I030 = I020 + IRPDPD(IRPBC,13)
C
C     I0T1 --- T1(B,N)
C     I000 --- W(NcdK) ordered (cN,dK)
C     I010 --- Z(B,c)
C     I020 --- Current Hbar element, read from disk, dist. size is B,c
C
      CALL ZERO(CORE(I010),IRPDPD(IRPBC,13))
      CALL ZERO(CORE(I020),IRPDPD(IRPBC,13))
      CALL ZERO(CORE(I030),IRPDPD(IRPBC,13))
C
      DO 260  IRPN=1,NIRREP
      IRPB  = DIRPRD(IRREPT,IRPN)
      IRPC  = DIRPRD(IRPB,IRPBC)
C
C     W(c,N) * T(B,N) (minus sign because of LWOFF+6).
C
      IOFFZ = I010 - 1 + ISYMOFF(IRPB,IRPBC,23)
C
      CALL XGEMM('N','T',VRT(IRPC,2),VRT(IRPB,1),
     &                               POP(IRPN,1),ONEM,
     &           CORE(IOFFW),VRT(IRPC,2),
     &           CORE(IOFFT),VRT(IRPB,1),ZILCH,
     &           CORE(IOFFZ),VRT(IRPC,2))
C
      IOFFT = IOFFT + VRT(IRPB,1) * POP(IRPN,1)
      IOFFW = IOFFW + VRT(IRPC,2) * POP(IRPN,1)
  260 CONTINUE
C
      CALL SYMTRA2(IRPBC,VRT(1,2),VRT(1,1),IRPDPD(IRPBC,13),1,
     &             CORE(I010),CORE(I030))
      KD = ISYMOFF(IRPD,IRPDK,18) + K - 1 + (D-1)*POP(IRPK,1)
      CALL GETLST(CORE(I020),KD,1,2,IRPDK,LISTZ)
      CALL  SAXPY(IRPDPD(IRPBC,13),ONE,CORE(I030),1,CORE(I020),1)
      CALL PUTLST(CORE(I020),KD,1,2,IRPDK,LISTZ)
C
  270 CONTINUE
  280 CONTINUE
  290 CONTINUE
  300 CONTINUE
      RETURN
      END
