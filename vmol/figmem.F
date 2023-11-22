      SUBROUTINE FIGMEM(MAXSHL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER XNEQV,AND,OR,EOR,DSTRT
      LOGICAL IFD1,IFD2,IFD3
C-----------------------------------------------------------------------
C     Parameters
C-----------------------------------------------------------------------
#include "baslims.par"
#include "mxatms.par"
      parameter (mxp2=maxprim*maxprim)
      parameter (khm=(nht*(nht+1))/2,nhl=mxp2*khm*khm)
      parameter (kwd=(nht+1)*(nht+2)*(nht+3)/6,nh4=4*nht-3)
C-----------------------------------------------------------------------
      DIMENSION LEOR8(64)
      COMMON /FLAGS/ IFLAGS(100)
      common /machsp/ iintln,ifltln,iintfp,ialone,ibitwd
      common /meminf/ maamraf0,memneed0,itotsiz,ihwm1,ihwm2
c      COMMON /ITY/ ITYA(28),ITYB(28),ITYC(28),ITYD(28),
c     1 KPQR,NSTA,NSTB,NSTC,NSTD,MUT,NUMSM(4),NBLOCK(4),
c     1 JA(8),JB(8),JC(8),JD(8),BUF(2400),IBUF(2400),NUTTY(4),N2TAPE(4),
c     X LEOR8(64),IFD3,MULA,MULB,MULC,MULD
      COMMON /INDX/ PC(512),DSTRT(8,MXCBF),NTAP,LU2,NRSS,NUCZ,ITAG,
     1 MAXLOP,MAXLOT,KMAX,NMAX,KHKT(7),MULT(8),ISYTYP(3),ITYPE(7,28),
     2 LAND(64),XNEQV(8,8),EOR(8,8),NPARSU(8),NPAR(8),MULNUC(Mxatms),
     3 NHKT(MXTNSH),MUL(MXTNSH),NUCO(MXTNSH),NRCO(MXTNSH),JSTRT(MXTNSH),
     4 NSTRT(MXTNSH),MST(MXTNSH),JRS(MXTNSH)
      INEQV4(I,J,K,L) = NEQV(NEQV(I,J) , NEQV(K,L))
c
      mim=maxlop-1
c
      DO 177 I=0,7
      DO 177 J=0,7
  177 LEOR8(8*I+J+1) = NEQV(I,J) * 8
c
c
      imemsz1=0
      ihwm1=0
      maamraf0=0
      memneed0=0
C
C loop over shells
C
      do 1 ishl1=1,maxshl
       mula=mul(ishl1)
      do 2 ishl2=1,ishl1
       mulb=mul(ishl2)
      do 3 ishl3=1,ishl1
       mulc=mul(ishl3)
       itop=ishl3
       if(ishl3.eq.ishl1)itop=ishl2
      do 4 ishl4=1,itop
       muld=mul(ishl4)

      KPQ=AND(MULA-1,MULB-1)
      MPQ=OR(MULA-1,MULB-1)
      KPQR=AND(KPQ,MULC-1)
      MPQR=OR(KPQ,MULC-1)
      MPQRS=OR(KPQR,MULD-1)
      MX = AND(KPQR,MULD-1)
      KPQRS = MX+1
      NMSS = MULT(MULA)*MULT(MULB)*MULT(MULC)*MULT(MULD)
      NMAS = NMSS/MULT(KPQRS)
       ifd1=ishl1.gt.ishl2
       ifd2=ishl3.gt.ishl4
       ifd3=ishl1.gt.ishl3.or.ishl2.ne.ishl4
       iang1=nhkt(ishl1)
       iang2=nhkt(ishl2)
       iang3=nhkt(ishl3)
       iang4=nhkt(ishl4)
       nang1=khkt(iang1)
       nang2=khkt(iang2)
       nang3=khkt(iang3)
       nang4=khkt(iang4)
       nprm1=nuco(ishl1) 
       nprm2=nuco(ishl2) 
       nprm3=nuco(ishl3) 
       nprm4=nuco(ishl4) 
       ncnt1=nrco(ishl1) 
       ncnt2=nrco(ishl2) 
       ncnt3=nrco(ishl3) 
       ncnt4=nrco(ishl4) 
       naqrs=nprm1*nprm2*nprm3*nprm4
       maa =nang1*nang2*nang3*nang4
       mraf=ncnt1*ncnt2*ncnt3*ncnt4
       mconlor=(iang1+iang2+iang3+iang4)*naqrs
       maamraf=maa*mraf
       nnn=iang1+iang2+iang3+iang4-3
       nnb=iang1*iang2*iang3*iang4*(nnn+1)/2
       npw=naqrs*max(2,nnn-1)
       nbb=nnb*naqrs
c
c now get lia and lina value
c
       lia=0
       lina=0
       do 100 indexi=1,khkt(iang1)
       ita=(itype(iang1,indexi)-1)*8
       do 110 indexj=1,khkt(iang2)
       if(.not.ifd1.and.indexj.gt.indexi)goto 110
       itb=(itype(iang2,indexj)-1)*8
       do 120 indexk=1,khkt(iang3)
       if(.not.ifd3.and.indexk.gt.indexi)goto 120
       itc=(itype(iang3,indexk)-1)*8
       do 130 indexl=1,khkt(iang4)
       if(.not.ifd2.and.indexl.gt.indexk)goto 130
       itd=(itype(iang4,indexl)-1)*8
      IF(AND(INEQV4(ITA,ITB,ITC,ITD)/8,MX) .NE. 0) LINA=LINA+NMAS
        lia=lia+1
       do 140 la=0,mim
        if(land(leor8(ita+la+1)+mula).ne.1)goto 140
       do 150 my=0,mim
        if(land(leor8(itb+my+1)+mulb).ne.1)goto 150
       do 160 ny=0,mim
        if(land(leor8(itc+ny+1)+mulc).ne.1)goto 160
        LMNY = NEQV(NEQV(LA,MY),NY)
        IF(LAND(LEOR8(ITD+LMNY+1)+MULD) .NE. 1)goto 170
        LINA=LINA+1
170    continue
160    continue
150    continue
140    continue
130    continue
120    continue
110    continue
100    continue
c       maxlina=64*lia
       maxlina=lina
c
c a 1 area
c
       isize1=maxlina*(mraf+1)+2*naqrs
       isize2a=nnn*naqrs+3*nbb+3*npw
       isize2b=1000000
       ihwm1  =max(isize1,ihwm1)
       imemsz1=max(isize1+max(isize2a,isize2b),imemsz1)
c
c a 2 area
c
       ifact=1
       if(iflags(62).eq.1)ifact=2
       memneed=max(4*naqrs,lia+3*naqrs+
     &             ifact*maamraf+max(mconlor,maamraf))+1
       maamraf0=max(maamraf0,maamraf)
       memneed0=max(memneed,memneed0)
       ihwm2=memneed0
c
4     continue
3     continue
2     continue
1     continue
c
CJDW 1/27/97. Work out the memory needed for ONELH.
C
      CALL CLCNRXLFT(NROWXLFT)
      IG = 40000
      IMEMONELH = 50 + 3*KWD + 3 + 3 + 3 + MXP2*8 + MXP2 + NHL +
     &             2 + 2 + 2 + IG + 100 + IG + 600 + 600 + IG + IG +
     &             NROWXLFT*8*3 + 3 + 14 + 3 + 11 + 7 + 3 + 51 + 3 +
     &             39 + 39 + IG + IG + 600
C
      ihwm1=(ihwm1+50)*iintfp
      ihwm2=(ihwm2+50)*iintfp
      imemsz1=max(imemsz1,imemonelh)*iintfp
      itotsiz=imemsz1+ihwm2 
c
      write(6,*)' required memory for a1 array ',imemsz1,' words '
      write(6,*)' required memory for a2 array ',ihwm2,' words '
      return
      end
