      SUBROUTINE CONLOX (F,SABB,SCDD,NPL,AINT,XLFT,XRGT,E1,
     $ SAVROT,TMUL,ROOT,QVV,NND,NNE,MAXLOT,I2,A2,I1,A1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DSTRT,AND,OR,EOR
C-----------------------------------------------------------------------
C     Parameters
C-----------------------------------------------------------------------
#include "baslims.par"
#include "mxatms.par"
      parameter (nh4=4*nht-3, nh2=nht+nht+1)
      parameter (mxp2=maxprim*maxprim)
      parameter (mxp4=mxp2*mxp2)
C-----------------------------------------------------------------------
      DIMENSION I2(20),A2(2),I1(50),A1(1)
      COMMON /ICOMP/ IBE(mxp4),K,NN
      common /icomp2/IBEM
      COMMON /REP/ NEWIND(MXCBF) , MSTOLD(8)
      common /meminf/ ijunk(3),ihwm1,ihwm2
      COMMON /MMMM/ MDUM(7),MRAF,JJUNK(4)
      COMMON /XA/ XAND(3,8,nh4)
      LOGICAL NPL(3,8,2),IFD1X,IFD2X,IFPL,IFLAG
      COMMON /FLAG/ IFLAG
      COMMON /ONE/ 
     1 ALP(mxp2),BET(mxp2),SS1(mxp2),CSS1(mxp2),PV(mxp2,8,3),CAB(mxp2),
     2 GAM(mxp2),DEL(mxp2),SS2(mxp2),CSS2(mxp2),QV(mxp2,8,3),CCD(mxp2)
      COMMON /TAB/ NNN
      COMMON/NULL_COM/FACT(nh4),RFACT(nh4),FACTM(nh4),RFACTM(nh4)
     & ,MAA,IFD1,IFD2,KCD,KBCD,NHCD,NHBCD,NNC,IFPL(3),
     & NN1,NN2,NHKTA,NHKTB,NHKTC,NHKTD,KHKTA,KHKTB,KHKTC,KHKTD,NNB
      COMMON /CON/ RAX,RAY,RAZ,RPRX,RPRY,RPRZ,DISTPQ,DISTRS,TOLR,TP52,
     &             ZJUNK(MXP2,16)
      DIMENSION SABB(2),SCDD(2),F(2),AINT(2),XLFT(NND,MAXLOT,3),
     $ XRGT(NNE,MAXLOT,3),E1(2),SAVROT(2),TMUL(2),ROOT(2),QVV(mxp2,8,3)
      COMMON /CONI/ NUCT(4),NRCT(4),JSTT(4),
     1 IFD1X,IFD2X,NUCAQ,NUCRS,NAQRS,NCQRS,KQ,KR,KS,IVA
      COMMON /INDX/ PC(512),DSTRT(8,MXCBF),NTAP,LU2,NRSS,NUCZ,ITAG,
     1 MAXLOP,MAXXXX,KMAX,NMAX,KHKT(7),MULT(8),ISYTYP(3),ITYPE(7,28),
     2 AND(8,8),OR(8,8),EOR(8,8),NPARSU(8),NPAR(8),MULNUC(Mxatms),
     3 NHKT(MXTNSH),MUL(MXTNSH),NUCO(MXTNSH),NRCO(MXTNSH),JSTRT(MXTNSH),
     4 NSTRT(MXTNSH),MST(MXTNSH),JRS(MXTNSH)
      COMMON /DAT/  ALPHA(MXTNPR),CONT(MXTNCC),CENT(3,MXTNSH),
     1              CORD(Mxatms,3),CHARGE(Mxatms),FMULT(8),TLA, TLC

#ifdef _TRACE_VMOL
c      write(*,*) '@VMOL: Entered file ',__FILE__
#endif

C     ROUTINE FOR GENERAL CONTRACTION CASE
      CALL MVPROD (E1,SABB,NUCAQ,SCDD,NUCRS)
      DO 70 I=1,NAQRS
      E1(I) = E1(I)*SAVROT(I)
   70 CONTINUE
      CALL CMPIN(E1,ROOT,TLC,NAQRS,NCQRS,F(2))
      IF(NCQRS .EQ. 0) RETURN
      NNN3 = MAX0(2,NNN-1)
      NBB=NNB*NCQRS
      NPW = NCQRS * MAX0(2,NNN-1)
      I1(8)=I1(7)+NNN*NAQRS
      I1(9)=I1(8)+3*NBB
      I1(10)=I1(9)
      I1(11)=I1(10)+3*NPW
C JFS
      I2(7)=I2(6)+(NHKTA+NHKTB+NHKTC+NHKTD)*NAQRS
      if(i2(7).gt.ihwm2)then
        write(*,*) '@CONLOX: VMOL memory error.'
        call errex
       endif
C
C peter szalay's contribution!
c
      call zero(a1(i1(10)),3*npw)
C
      IF(I1(11).GE.I1(50).OR.I2(7).GT.I2(20)) CALL LARM(I1,I2)
      DO 40 I=1,3
      IF(IFPL(I)) GO TO 40
      CALL PQ (E1,PV(1,KQ,I),NUCAQ,QVV(1,IVA,I),NUCRS,XAND(I,KR,2))
      CALL CMPR1 (E1,A1(I1(10)+(I-1)*NPW))
c     call gather (ibem,A1(I1(10)+(I-1)*NPW),E1,ibe)
   40 CONTINUE
      CALL CMPR1(TMUL,E1)
c     call gather(ibem,E1,TMUL,IBE)
      IGAM=MAX0(5,NNN)*NCQRS
      CALL TGCAL (A1(I1(10)),F,NCQRS,NNN3,IFPL)
      DO 130 LQS=1,NCQRS
  130 F(IGAM+LQS) = F(LQS)*E1(LQS)
cold  CALL VGAMMA (F,F,NNN-1,NCQRS)
      CALL VGAMMA2(F,F(1+IGAM),F(1+IGAM+NCQRS),F(1+IGAM+2*NCQRS),
     &             NNN-1,NCQRS)
      DO 140 LQS=1,NCQRS
  140 F(LQS) = F(LQS)*ROOT(LQS)
      IB = 0
      DO 160 I=2,NNN
      IB = IB + NCQRS
      DO 150 LQS=1,NCQRS
      ROOT(LQS) = -4.*ROOT(LQS)*E1(LQS)
      F(LQS+IB) = F(LQS+IB)*ROOT(LQS)
  150 CONTINUE
  160 CONTINUE
      CALL ZERO (A1(I1(8)),3*NBB)
      CALL DORBAX (F,NPL(1,KQ,1),NPL(1,IVA,2),A1(I1(8)),A1(I1(10)),NBB,
     & NPW,XLFT,XRGT,A2(I2(6)),NAQRS,NND,NNE,MAXLOT)
      CALL DORBBR (A1(I1(7)),A1(I1(9)),A1(I1(8)),NBB,AINT)
#ifdef _TRACE_VMOL
c      write(*,*) '@VMOL: Leaving file ',__FILE__
#endif
      RETURN
      END
