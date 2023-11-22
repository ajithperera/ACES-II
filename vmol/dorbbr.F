      SUBROUTINE DORBBR (HLP,HLP2,BBXS,NBB,AINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DSTRT,AND,OR,EOR
C-----------------------------------------------------------------------
C     Parameters
C-----------------------------------------------------------------------
#include "baslims.par"
#include "mxatms.par"
      parameter (nh4=4*nht-3, nh2=nht+nht+1)
      parameter (nhq=nht*nht*nht*nht)
      parameter (mxp2=maxprim*maxprim)
C-----------------------------------------------------------------------
      DIMENSION BBXS(NBB,3),HLP(2),HLP2(2),AINT(2)
      LOGICAL  LG1,LG2,LG3,IFD1,IFD2,IFD3,IFD1X,IFD2X,IFLAG,t3,lharm
      COMMON /FLAG/ IFLAG
      COMMON /MMMM/ MDUM(7),MRAF,JJUNK(4)
      COMMON /INDX/ PC(512),DSTRT(8,MXCBF),NTAP,LU2,NRSS,NUCZ,ITAG,
     1 MAXLOP,MAXXXX,KMAX,NMAX,KHKT(7),MULT(8),ISYTYP(3),ITYPE(7,28),
     2 IVAND(8,8),IVOR(8,8),IVEOR(8,8),NPARSU(8),NPAR(8),MULNUC(Mxatms),
     3 NHKT(MXTNSH),MUL(MXTNSH),NUCO(MXTNSH),NRCO(MXTNSH),JSTRT(MXTNSH),
     4 NSTRT(MXTNSH),MST(MXTNSH),JRS(MXTNSH)
c      COMMON /INDX/ PC(512),DSTRT(8,500),NTAP,LU2,NRSS,NUCZ,ITAG
      COMMON /ITY/ ITYA(28),ITYB(28),ITYC(28),ITYD(28),
     1 KPQR,NSTA,NSTB,NSTC,NSTD,MUT,NUMSM(4),NBLOCK(4),
     1 JA(8),JB(8),JC(8),JD(8),BUF(2400),IBUF(2400),NUTTY(4),N2TAPE(4),
     X LEOR8(64),IFD3,MULA,MULB,MULC,MULD
      COMMON /TAB/ NNN
      COMMON/NULL_COM/FACT(nh4),RFACT(nh4),FACTM(nh4),RFACTM(nh4)
     & ,MAA,IFD1,IFD2,KCD,KBCD,NHCD,NHBCD,NNC,LG1,LG2,LG3,
     & NN1,NN2,NHKTA,NHKTB,NHKTC,NHKTD,KHKTA,KHKTB,KHKTC,KHKTD,NJUNK
      COMMON /DAT/  ALPHA(MXTNPR),CONT(MXTNCC),CENT(3,MXTNSH),
     1              CORD(Mxatms,3),CHARGE(Mxatms),FMULT(8),TLA, TLC
      COMMON /CONI/ NUCT(4),NRCT(4),JSST(4),
     1 IFD1X,IFD2X,NUCAQ,NUCRS,NAQRS,NCQRS,KQ,KR,KS,IVA
      COMMON /DOR/ INDXY(nhq),INDZ(nhq)
      COMMON /HARMON/ LHARM

#ifdef _TRACE_VMOL
c      write(*,*) '@VMOL: Entered file ',__FILE__
#endif

      NUM=0
      numa=0
      DO 500 L1=NHKTA-1,0,-1
      DO 500 N1=0,NHKTA-1-L1
      numa=numa+1
      M1 = NHKTA-1-L1-N1
      numb=0
      DO 400 L2=NHKTB-1,0,-1
      DO 400 N2=0,NHKTB-1-L2
      numb=numb+1
      if(numb.gt.numa.and..not.ifd1.and..not.lharm) GO TO 400
      M2 = NHKTB-1-L2-N2
      LL12=L1+L2+1
      MM12=M1+M2+1
      IX2=L1*NHBCD+L2*NHCD
      IY2=M1*NHBCD+M2*NHCD
      IZ2=N1*NHBCD+N2*NHCD
      numc=0
      DO 300 L3=NHKTC-1,0,-1
      DO 300 N3=0,NHKTC-1-L3
      numc=numc+1
      if(.not.ifd3 .and. numc.gt.numa .and. .not.lharm) go to 300
      t3=.not.ifd3 .and. numc.eq.numa
      M3 = NHKTC-1-L3-N3
      L123=LL12+L3
      M123=MM12+M3
      IX3=IX2+L3*NHKTD
      IY3=IY2+M3*NHKTD
      IZ3=IZ2+N3*NHKTD
      numd=0
      DO 200 L4=NHKTD-1,0,-1
      DO 200 N4=0,NHKTD-1-L4
      numd=numd+1
      if (.not.lharm .and. (numd.gt.numc  .and.  .not.ifd2)) go to 200
c     if (.not.lharm .and.  (numd.gt.numb .and. t3)) go to 200
      NUM = NUM + 1
      JT1 = L123 + L4
      M4 = NHKTD-1-L4-N4
      IF (LG1) THEN
      IF (AND(JT1,1) .EQ. 0) GO TO 200
      JT1H = JT1/2 + 1
      ELSE
      JT1H = JT1
      ENDIF
      JT2 = M123 + M4
      IF (LG2) THEN
      IF (AND(JT2,1) .EQ. 0) GO TO 200
      JT2H = JT2/2 + 1
      ELSE
      JT2H = JT2
      ENDIF
      IF(LG3 .AND. AND(1,NNN+2-JT1-JT2).EQ.0) GO TO 200
      CALL RCON (HLP2,BBXS(INDXY(IX3+L4+1),1),
     &                BBXS(INDXY(IY3+M4+1),2),
     &            BBXS(INDZ (IZ3+N4+1),3),
     &       HLP,NCQRS,JT1H,JT2H)
      CALL ZERO (HLP,NAQRS)
      CALL EXPND (HLP,HLP2)
      CALL TRANS (5-ITAG,CONT,JSST(ITAG),HLP,NUCT(ITAG),NRCT(ITAG),ISTO)
      IF(ITAG .EQ. 1) THEN
      CALL VM_MOVE(HLP(ISTO+1),AINT(MRAF*(NUM-1)+1),MRAF)
      ELSE IF (ITAG.EQ. 2) THEN
      IA=(NUCZ-1)/NUCT(2)+1
      NIN=MRAF/NRCT(1)
      DO 60 IR=0,NRCT(1)-1
      IE=IR*NIN + (NUM-1)*MRAF
      CO=CONT(JSST(1)+IR*NUCT(1)+IA)
      DO 60 I=1,NIN
   60 AINT(IE+I) = AINT(IE+I) + HLP(I+ISTO)*CO
      ELSE
      NIN=MRAF/(NRCT(1)*NRCT(2))
      IA=(NUCZ-1)/NUCT(2)
      IB=NUCZ-IA*NUCT(2)
      DO 70 IR=0,NRCT(1)-1
      CA=CONT(JSST(1)+IR*NUCT(1)+IA+1)
      DO 70 IS=0,NRCT(2)-1
      IE=(IR*NRCT(2)+IS)*NIN + (NUM-1)*MRAF
      CO=CA*CONT(JSST(2)+IS*NUCT(2)+IB)
      DO 70 I=1,NIN
   70 AINT(IE+I) = AINT(IE+I) + CO*HLP(ISTO+I)
      ENDIF
  200 CONTINUE
  300 CONTINUE
  400 CONTINUE
  500 CONTINUE

#ifdef _TRACE_VMOL
c      write(*,*) '@VMOL: Leaving file ',__FILE__
#endif

      RETURN
      END
