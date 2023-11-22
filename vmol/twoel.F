      SUBROUTINE TWOEL(ICORE,JCORE,NSYMOP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 TEXT(2),NTYPI(4)
      INTEGER AND,OR,EOR,DSTRT 
C-----------------------------------------------------------------------
C     Parameters
C-----------------------------------------------------------------------
#include "baslims.par"
#include "mxatms.par"
      parameter (nh4=4*nht-3)
      parameter (mxp2=maxprim*maxprim)
C-----------------------------------------------------------------------
      DIMENSION ICORE(*),JCORE(*)
      COMMON /XA/ XAND(3,8,nh4)
      COMMON /FLAG/ IFLAG
      LOGICAL IFLAG
      COMMON /MMMM/ M2,M3,LINA,MUM,MIM,NMAS,MX,MRAF,JJUNK(4)
      COMMON /TAB/ NNN
      COMMON /ONE/ 
     1   ALP(mxp2),BET(mxp2),S1(mxp2),CSS1(mxp2),PV(mxp2,8,3),CAB(mxp2),
     2   GAM(mxp2),DEL(mxp2),S2(mxp2),CSS2(mxp2),QV(mxp2,8,3),CCD(mxp2)
      LOGICAL  IFD1,IFD2,IFD3,IFPL,
     1 IFD1X,IFD2X,NPL(3,8,2)
      COMMON/NULL_COM/FACT(nh4),RFACT(nh4),FACTM(nh4),RFACTM(nh4)
     & ,MAA,IFD1,IFD2,KCD,KBCD,NHCD,NHBCD,NNC,IFPL(3),
     & NN1,NN2,NHKTA,NHKTB,NHKTC,NHKTD,KHKTA,KHKTB,KHKTC,KHKTD,NNB
      COMMON /CON/ RAX,RAY,RAZ,RPRX,RPRY,RPRZ,DISTPQ,DISTRS,TOLR,TP52,
     1 SABB(mxp2,8),SCDD(mxp2,8)
      COMMON /CONI/ NUCR,NUCS,NUCA,NUCQ,
     1 NRCR,NRCS,NRCA,NRCQ,  JSTR,JSTS,JSTA,JSTQ,
     1 IFD1X,IFD2X,NUCAQ,NUCRS,NAQRS,NCQRS,KQ,KR,KS,IVA
      COMMON /INDX/ PC(512),DSTRT(8,MXCBF),NTAP,LU2,NRSS,NUCZ,ITAG,
     1 MAXLOP,MAXLOT,KMAX,NMAX,KHKT(7),MULT(8),ISYTYP(3),ITYPE(7,28),
     2 AND(8,8),OR(8,8),EOR(8,8),NPARSU(8),NPAR(8),MULNUC(Mxatms),
     3 NHKT(MXTNSH),MUL(MXTNSH),NUCO(MXTNSH),NRCO(MXTNSH),JSTRT(MXTNSH),
     4 NSTRT(MXTNSH),MST(MXTNSH),JRS(MXTNSH)
      COMMON /DAT/  ALPHA(MXTNPR),CONT(MXTNCC),CENT(3,MXTNSH),
     1              CORD(Mxatms,3),CHARGE(Mxatms),FMULT(8),TLA, TLC
      COMMON /ITY/ ITYA(28),ITYB(28),ITYC(28),ITYD(28),
     1 KPQR,NSTA,NSTB,NSTC,NSTD,MUT,NUMSM(4),NBLOCK(4),
     1 JA(8),JB(8),JC(8),JD(8),BUF(2400),IBUF(2400),NUTTY(4),N2TAPE(4),
     X LEOR8(64),IFD3,MULA,MULB,MULC,MULD
      DATA TEXT/'********','TWOELSUP'/
      DATA NTYPI /' I I I I',' I J I J',' I I J J',' I J K L'/
      N2TAPE(1)=2
      N2TAPE(2)=81
      N2TAPE(3)=82
      N2TAPE(4)=83
      N2TAP = N2TAPE(1)
      JCORE(1)=1
      PI=2.0*DACOS(0.0D00)
      MIM = MAXLOP-1
      M2 = MAXLOP*MAXLOP
      M3 = M2*MAXLOP - 1
      PIROT=2.0/SQRT(PI)
C
C Originally the 8889 loop was executed 1 through 4 reagrdless
C of the status of the files. However, for non symmetry (C1) 
C cases, except IIII, the other files (IIJJ, IJIJ, IJKL) do 
C not get open in READIN.F. This caused problems on SGI
C machines (noted by Marshall Cory). In the following fix 
C we no longer try to write to files that are not opened.
C Ajith Perera 01/2002
C
      IF (NSYMOP .EQ. 0) THEN
          ITARGET = 1
      ELSE IF (NSYMOP .EQ. 1) THEN
          ITARGET = 3 
      ELSE 
          ITARGET = 4 
      ENDIF 
C
      DO 8889 I = 1, ITARGET
      NUTTY(I) = 0
      NUMSM(I) = 0
      NBLOCK(I) = 0
      N2TAPL = N2TAPE(I)
      WRITE(N2TAPL) TEXT(1) 
      WRITE(N2TAPL) TEXT(2) 
8889  CONTINUE
C
      TP52=PI*PI*PI*PIROT
      DO 177 I=0,7
      DO 177 J=0,7
  177 LEOR8(8*I+J+1) = NEQV(I,J) * 8
      NUT=0
      DO 1744 IA=1,KMAX
      JSTA=JRS(IA)
      NUCA=NUCO(IA)
      NRCA = NRCO(IA)
      MULA=MUL(IA)
      NSTA=NSTRT(IA)
      NHKTA=NHKT(IA)
      KHKTA=KHKT(NHKTA)
      DO 452 NA=1,KHKTA
  452 ITYA(NA)=(ITYPE(NHKTA,NA)-1)*8
      DO 1744 IB=1,IA
      MULB=MUL(IB)
      NHKTB=NHKT(IB)
      KHKTB=KHKT(NHKTB)
      NUCQ=NUCO(IB)
      NRCQ = NRCO(IB)
      JSTQ=JRS(IB)
      DO 453 NB=1,KHKTB
  453 ITYB(NB)=(ITYPE(NHKTB,NB)-1)*8
      IFD1=IA.GT.IB
      NSTB=NSTRT(IB)
      NUCAQ=NUCA*NUCQ
      NN1 = NHKTA+NHKTB-1
      NN11=NUCAQ*NHKTA*NHKTB*(NN1+1)/2
      JCORE(2)=JCORE(1)+3*MAXLOT*NN11
      CALL EXPA(IA,IB,SABB,FACT,JCORE(50+JCORE(1)),ALP,BET,S1,
     & CSS1,PV,NPL,CAB,
     $ NN11,MAXLOT)
      DO 8 IC=1,IA
      NHKTC=NHKT(IC)
      NUCR=NUCO(IC)
      NRCR=NRCO(IC)
      JSTR=JRS(IC)
      MULC=MUL(IC)
      KHKTC=KHKT(NHKTC)
      DO 454 NC=1,KHKTC
  454 ITYC(NC)=(ITYPE(NHKTC,NC)-1)*8
      NSTC=NSTRT(IC)
      IDS=IC
      IF(IC.EQ.IA) IDS = IB
      DO 8 ID=1,IDS
      MULD=MUL(ID)
      NHKTD=NHKT(ID)
      NUCS=NUCO(ID)
      NRCS=NRCO(ID)
      JSTS=JRS(ID)
      KHKTD=KHKT(NHKTD)
      DO 455 ND=1,KHKTD
  455 ITYD(ND)=(ITYPE(NHKTD,ND)-1)*8
      IFD2=IC.GT.ID
      IFD3=IA.GT.IC.OR.IB.NE.ID
      NN2=NHKTC+NHKTD-1
      KCD=KHKTC*KHKTD
      KBCD=KCD*KHKTB
      NHCD=NHKTC*NHKTD
      NHBCD=NHKTB*NHCD
      NNC=NN2+NHKTB-1
      NSTD=NSTRT(ID)
      MAA=KHKTA*KHKTB*KHKTC*KHKTD
      MRAF = NRCA*NRCQ*NRCR*NRCS
      NNN=NN1+NN2-1
      NN22=NUCR*NUCS*NHCD*(NN2+1)/2
      JCORE(3)=JCORE(2)+3*MAXLOT*NN22
      CALL DRIVE(ICORE(1),ICORE(21),JCORE(1),JCORE(51),IA,IC,ID,NPL)
    8 CONTINUE
 1744 CONTINUE
C...  CLOSING SEQUENCE
      IBOFF = 0
      NUMSMT = 0
      NBLOCQ = 0
      DO 8898 I = 1, ITARGET
      N2TAP = N2TAPE(I)
      IF (NUTTY(I) .NE. 0) THEN
        NBLOCK(I) = NBLOCK(I) + 1
        call wrseq(buf(1+iboff),ibuf(1+iboff),nutty(i),600,n2tap)
c        WRITE(N2TAP) (BUF(IBOFF+J),J = 1,600),
c     1               (IBUF(IBOFF+J),J = 1,600), NUTTY(I)
      ENDIF
C...  WRITE EMPTY BLOCK TO INDICATE END
      NUT=-1
      call wrseq(buf(1+iboff),ibuf(1+iboff),nut,600,n2tap)
c      WRITE(N2TAP) (BUF(IBOFF+J),J = 1,600),
c     1             (IBUF(IBOFF+J),J = 1,600), NUT
      END FILE N2TAP
      REWIND N2TAP
      IF (NUMSM(I) .NE. 0)THEN
       WRITE(6,2000) NUMSM(I), NTYPI(I)
2000   FORMAT(T3,'@TWOEL-I, ',I10,' integrals of symmetry type ',A)
      ENDIF
      IBOFF = IBOFF + 600
      NUMSMT = NUMSMT + NUMSM(I)
      NBLOCQ = NBLOCQ + NBLOCK(I)
8898  CONTINUE
      WRITE(6,2100) NUMSMT
2100  FORMAT(T3,'@TWOEL-I, Total number of 2-e integrals ',I10,'.')
      RETURN
      END
