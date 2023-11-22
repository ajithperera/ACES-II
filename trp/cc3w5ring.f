      SUBROUTINE CC3W5RING(CORE,MAXCOR,IUHF)
C
C This subroutine is used to compute the contribution
C
C     Z(bc,dk) = - \sum_n [ t(c,n) * <nb||kd> + t(b,n) * <nc||kd> ]
C     Z        =            T * W
C
C to the Hbar element. It is called in CC3 calculations. A similar
C contraction is performed by W5RING, but CC3W5RING needs only integrals.
C Once CC3W5RING is debugged, it may well be replaced by a generalized
C W5RING.
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
      LISTT1 =  90
      I0T1(1) = 1
      I0T1(2) = I0T1(1) + NT(1)
      CALL GETLST(CORE(I0T1(1)),1,1,1,1     ,LISTT1)
      CALL GETLST(CORE(I0T1(2)),1,1,1,1+IUHF,LISTT1)
C
C     Z(Bc,Dk) = -t(c,n) * <nB||kD> + t(B,N) * <Nc||kD>
C
      LISTZ  = 130
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
      IF(NEED.GT.MAXCOR) CALL INSMEM('CC3W5RING',NEED,MAXCOR)
C
C --- Pick up <nB||kD> ordered (Dn,Bk). Reorder to (Bn,Dk) ---
C
      CALL GETALL(CORE(I010),ISIZW,1,LISTW)
      CALL SSTGEN(CORE(I010),CORE(I000),ISIZW,VRT(1,1),POP(1,2),
     &            VRT(1,1),POP(1,2),CORE(I020),1,'3214')
C
      IOFFW = I000
C
      DO  50 IRPDK=1,NIRREP
C
      IRPBC = IRPDK
C
      DK = 1
C
      DO  40  IRPK=1,NIRREP
      IRPD  = DIRPRD(IRPDK,IRPK)
C
      DO  30     K=1,POP(IRPK,2)
      DO  20     D=1,VRT(IRPD,1)
C
      IOFFT = I0T1(2)
      IOFFZ = I010
C
      I020 = I010 + IRPDPD(IRPBC,13)
      I030 = I020 + IRPDPD(IRPBC,13)
C
C     I0T1 --- T1(c,n)
C     I000 --- <nB||kD> ordered (bn,dk)
C     I010 --- Z(B,c)
C     I020 --- Current Hbar element, read from disk, dist. size is B,c
C
      DO  10  IRPN=1,NIRREP
      IRPC  = IRPN
      IRPB  = DIRPRD(IRPC,IRPBC)
C
C     -W(B,n) * T(c,n)
C
      CALL XGEMM('N','T',VRT(IRPB,1),VRT(IRPC,2),
     &                               POP(IRPN,2),ONEM,
     &           CORE(IOFFW),VRT(IRPB,1),
     &           CORE(IOFFT),VRT(IRPC,2),ZILCH,
     &           CORE(IOFFZ),VRT(IRPB,1))
C
      IOFFT = IOFFT + VRT(IRPC,2) * POP(IRPN,2)
      IOFFW = IOFFW + VRT(IRPB,1) * POP(IRPN,2)
      IOFFZ = IOFFZ + VRT(IRPB,1) * VRT(IRPC,2)
   10 CONTINUE
C
      CALL GETLST(CORE(I020),DK,1,2,IRPDK,LISTZ)
      CALL  SAXPY(IRPDPD(IRPBC,13),ONE,CORE(I010),1,CORE(I020),1)
      CALL PUTLST(CORE(I020),DK,1,2,IRPDK,LISTZ)
C
      DK = DK + 1
   20 CONTINUE
   30 CONTINUE
   40 CONTINUE
   50 CONTINUE
C
C     Z(Bc,Dk) = t(B,N) * <Nc||kD>
C
      LISTZ  = 130
C
      LISTW  =  21 + IUHF
      ISIZW = ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
C
      I000    = I0T1(2) + NT(2)
      I010    = I000    + ISIZW
      I020    = I010    + ISIZW
C     Note I020 length can be cut.
      I030    = I020    + ISIZW
      NEED    = I030 * IINTFP
      IF(NEED.GT.MAXCOR) CALL INSMEM('CC3W5RING',NEED,MAXCOR)
C
C --- Pick up <Nc||kD> = -(ND|ck) = -<Dc|Nk>. In UHF integral is on list
C --- 22 and is ordered (cN,Dk). In RHF <Dc|Nk>=<cD|kN>=<Cd|Kn> and we
C --- can get the integral from list 21 with (Cn,dK)=(cN,Dk) order.
C
      CALL GETALL(CORE(I000),ISIZW,1,LISTW)
C
      IOFFW = I000
C
      DO 100 IRPDK=1,NIRREP
C
      IRPBC = IRPDK
C
      DK = 1
C
      DO  90  IRPK=1,NIRREP
      IRPD  = DIRPRD(IRPDK,IRPK)
C
      DO  80     K=1,POP(IRPK,2)
      DO  70     D=1,VRT(IRPD,1)
C
      IOFFT = I0T1(1)
      IOFFZ = I010
C
      I020 = I010 + IRPDPD(IRPBC,13)
      I030 = I020 + IRPDPD(IRPBC,13)
C
C     I0T1 --- T1(B,N)
C     I000 --- <Nc||kD> ordered (cN,Dk)
C     I010 --- Z(B,c)
C     I020 --- Current Hbar element, read from disk, dist. size is B,c
C
      DO  60  IRPN=1,NIRREP
      IRPB  = IRPN
      IRPC  = DIRPRD(IRPB,IRPBC)
C
C     +W(c,N) * T(B,N)
C
      CALL XGEMM('N','T',VRT(IRPC,2),VRT(IRPB,1),
     &                               POP(IRPN,1),ONEM,
     &           CORE(IOFFW),VRT(IRPC,2),
     &           CORE(IOFFT),VRT(IRPB,1),ZILCH,
     &           CORE(IOFFZ),VRT(IRPC,2))
C
      IOFFT = IOFFT + VRT(IRPB,1) * POP(IRPN,1)
      IOFFW = IOFFW + VRT(IRPC,1) * POP(IRPN,1)
      IOFFZ = IOFFZ + VRT(IRPB,1) * VRT(IRPC,2)
   60 CONTINUE
C
      CALL SYMTRA2(IRPBC,VRT(1,2),VRT(1,1),IRPDPD(IRPBC,13),1,
     &             CORE(I010),CORE(I030))
      CALL GETLST(CORE(I020),DK,1,2,IRPDK,LISTZ)
      CALL  SAXPY(IRPDPD(IRPBC,13),ONE,CORE(I030),1,CORE(I020),1)
      CALL PUTLST(CORE(I020),DK,1,2,IRPDK,LISTZ)
C
      DK = DK + 1
   70 CONTINUE
   80 CONTINUE
   90 CONTINUE
  100 CONTINUE
C
      IF(IUHF.EQ.0) RETURN
C
      DO 200 ISPIN=1,2
C
      LISTZ = 126 + ISPIN
      LISTW =  22 + ISPIN
      LISTT1 = 90
C
      ISIZW = ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
C
      I000    = I0T1(2) + NT(2)
      I010    = I000    + ISIZW
      I020    = I010    + ISIZW
C     Note I020 length can be cut.
      I030    = I020    + ISIZW
      NEED    = I030 * IINTFP
      IF(NEED.GT.MAXCOR) CALL INSMEM('CC3W5RING',NEED,MAXCOR)
C
C --- Pick up <nb||kd> ordered (dn,bk). Reorder to (bn,dk) ---
C
      CALL GETALL(CORE(I010),ISIZW,1,LISTW)
      CALL SSTGEN(CORE(I010),CORE(I000),ISIZW,VRT(1,ISPIN),POP(1,ISPIN),
     &            VRT(1,ISPIN),POP(1,ISPIN),CORE(I020),1,'3214')
C
      IOFFW = I000
C
      DO 150 IRPDK=1,NIRREP
C
      IRPBC = IRPDK
C
      DK = 1
C
      DO 140  IRPK=1,NIRREP
      IRPD  = DIRPRD(IRPDK,IRPK)

      DO 130     K=1,POP(IRPK,ISPIN)
      DO 120     D=1,VRT(IRPD,ISPIN)
C
      IOFFT = I0T1(ISPIN)
      IOFFZ = I010
C
      I020 = I010 + IRPDPD(IRPBC,18+ISPIN)
      I030 = I020 + IRPDPD(IRPBC,   ISPIN)
C
C     I0T1 --- T1(c,n)
C     I000 --- <nb||kd> ordered (bn,dk)
C     I010 --- Z(b,c); after ASSYM2 becomes Z(b<c)
C     I020 --- Current Hbar element, read from disk, dist. size is b<c
C
      DO 110  IRPN=1,NIRREP
      IRPC  = IRPN
      IRPB  = DIRPRD(IRPC,IRPBC)
C
C     -W(b,n) * T(c,n)
C
      CALL XGEMM('N','T',VRT(IRPB,ISPIN),VRT(IRPC,ISPIN),
     &                                   POP(IRPN,ISPIN),ONEM,
     &           CORE(IOFFW),VRT(IRPB,ISPIN),
     &           CORE(IOFFT),VRT(IRPC,ISPIN),ZILCH,
     &           CORE(IOFFZ),VRT(IRPB,ISPIN))
C
      IOFFT = IOFFT + VRT(IRPC,ISPIN) * POP(IRPN,ISPIN)
      IOFFW = IOFFW + VRT(IRPB,ISPIN) * POP(IRPN,ISPIN)
      IOFFZ = IOFFZ + VRT(IRPB,ISPIN) * VRT(IRPC,ISPIN)
  110 CONTINUE
C
      CALL ASSYM2(IRPBC,VRT(1,ISPIN),1,CORE(I010))
      CALL GETLST(CORE(I020),DK,1,2,IRPDK,LISTZ)
      CALL  SAXPY(IRPDPD(IRPBC,ISPIN),ONE,CORE(I010),1,CORE(I020),1)
      CALL PUTLST(CORE(I020),DK,1,2,IRPDK,LISTZ)
C
      DK = DK + 1
  120 CONTINUE
  130 CONTINUE
  140 CONTINUE
  150 CONTINUE
  200 CONTINUE
C
C     Z(Bc,Kd) = -t(c,n) * <nB||Kd> + t(B,N) * <Nc||Kd>
C
      LISTZ  = 129
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
      IF(NEED.GT.MAXCOR) CALL INSMEM('CC3W5RING',NEED,MAXCOR)
C
C --- Pick up <nB||Kd> ordered (Bn,dK). ---
C
      CALL GETALL(CORE(I000),ISIZW,1,LISTW)
C
      IOFFW = I000
C
      DO 250 IRPDK=1,NIRREP
C
      IRPBC = IRPDK
C
      DO 240  IRPK=1,NIRREP
      IRPD  = DIRPRD(IRPDK,IRPK)
C
      DO 230     K=1,POP(IRPK,1)
      DO 220     D=1,VRT(IRPD,2)
C
      IOFFT = I0T1(2)
      IOFFZ = I010
C
      I020 = I010 + IRPDPD(IRPBC,13)
      I030 = I020 + IRPDPD(IRPBC,13)
C
C     I0T1 --- T1(c,n)
C     I000 --- <nB||Kd> ordered (bn,dk)
C     I010 --- Z(B,c)
C     I020 --- Current Hbar element, read from disk, dist. size is B,c
C
      DO 210  IRPN=1,NIRREP
      IRPC  = IRPN
      IRPB  = DIRPRD(IRPC,IRPBC)
C
C     -W(B,n) * T(c,n)
C
      CALL XGEMM('N','T',VRT(IRPB,1),VRT(IRPC,2),
     &                               POP(IRPN,2),ONEM,
     &           CORE(IOFFW),VRT(IRPB,1),
     &           CORE(IOFFT),VRT(IRPC,2),ZILCH,
     &           CORE(IOFFZ),VRT(IRPB,1))
C
      IOFFT = IOFFT + VRT(IRPC,2) * POP(IRPN,2)
      IOFFW = IOFFW + VRT(IRPB,1) * POP(IRPN,2)
      IOFFZ = IOFFZ + VRT(IRPB,1) * VRT(IRPC,2)
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
C     Z(Bc,Kd) = t(B,N) * <Nc||Kd>
C
      LISTZ  = 129
C
      LISTW  =  26
      ISIZW = ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
C
      I000    = I0T1(2) + NT(2)
      I010    = I000    + ISIZW
      I020    = I010    + ISIZW
C     Note I020 length can be cut.
      I030    = I020    + ISIZW
      NEED    = I030 * IINTFP
      IF(NEED.GT.MAXCOR) CALL INSMEM('CC3W5RING',NEED,MAXCOR)
C
C --- Pick up <Nc||Kd>, ordered (dN,cK). Reorder to (cN,dK).
C
      CALL GETALL(CORE(I010),ISIZW,1,LISTW)
      CALL SSTGEN(CORE(I010),CORE(I000),ISIZW,VRT(1,2),POP(1,1),
     &            VRT(1,2),POP(1,1),CORE(I030),1,'3214')
C
      IOFFW = I000
C
      DO 300 IRPDK=1,NIRREP
C
      IRPBC = IRPDK
C
      DO 290  IRPK=1,NIRREP
      IRPD  = DIRPRD(IRPDK,IRPK)
C
      DO 280     K=1,POP(IRPK,1)
      DO 270     D=1,VRT(IRPD,2)
C
      IOFFT = I0T1(1)
      IOFFZ = I010
C
      I020 = I010 + IRPDPD(IRPBC,13)
      I030 = I020 + IRPDPD(IRPBC,13)
C
C     I0T1 --- T1(B,N)
C     I000 --- <Nc||kD> ordered (cN,Dk)
C     I010 --- Z(B,c)
C     I020 --- Current Hbar element, read from disk, dist. size is B,c
C
      DO 260  IRPN=1,NIRREP
      IRPB  = IRPN
      IRPC  = DIRPRD(IRPB,IRPBC)
C
C     +W(c,N) * T(B,N)
C
      CALL XGEMM('N','T',VRT(IRPC,2),VRT(IRPB,1),
     &                               POP(IRPN,1),ONEM,
     &           CORE(IOFFW),VRT(IRPC,2),
     &           CORE(IOFFT),VRT(IRPB,1),ZILCH,
     &           CORE(IOFFZ),VRT(IRPC,2))
C
      IOFFT = IOFFT + VRT(IRPB,1) * POP(IRPN,1)
      IOFFW = IOFFW + VRT(IRPC,2) * POP(IRPN,1)
      IOFFZ = IOFFZ + VRT(IRPB,1) * VRT(IRPC,2)
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
