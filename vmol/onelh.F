      SUBROUTINE ONELH(I1,NXYZ,LYDI,LOVE,LMN,SABB,SAAB,HLP1,JJRS,NUC,
     &                 NRC,HVEKT,F,HINT,HBUF,JBUF,TSFSTA,TSAV,XLFT,
     &                 PP,TF,MNO,CP,XK,PYA,POWAB,MAKK,AIRU,SFX,TINT,
     &                 TVEKT,TBUF,ECPVECT,IG,NHL,NROWXLFT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,EOR,DSTRT
C-----------------------------------------------------------------------
C     Parameters
C-----------------------------------------------------------------------
#include "baslims.par"
#include "mxatms.par"
      PARAMETER (KWD=(NHT+1)*(NHT+2)*(NHT+3)/6)
      parameter (nh4=4*NHT-3)
      parameter (mxp2=maxprim*maxprim)
      parameter (khm=nht*(nht+1)/2)
C-----------------------------------------------------------------------
cch---------------------------------------------------------------------
c     ecp: array of ecpints
cch---------------------------------------------------------------------
c      common /ecpint/ ecpvekt(80000)
      dimension ecpvect(1)
      common /ecpsts/ iecp
cch---------------------------------------------------------------------
      COMMON /XA/ XAND(3,8,nh4)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /ONE/ 
     1 ALP(mxp2),BET(mxp2),S1(mxp2),CSS1(mxp2),PV(mxp2,8,3),CAB(mxp2),
     2 GAM(mxp2),DELT(mxp2),S2(mxp2),CSS2(mxp2),QV(mxp2,8,3),CCD(mxp2)
c      COMMON /ONE/ 
c     1   ALP(mxp2),BET(mxp2),S1(mxp2),CSS1(mxp2),PV(mxp2,8,3),CAB(mxp2)
      LOGICAL NPL(3,8), HARMT
      DIMENSION I1(50),NXYZ(KWD,3),LYDI(3),LOVE(3),LMN(3),
     $  SABB(MXP2,8), SAAB(mxp2),HLP1(NHL),
     & JJRS(2),NUC(2),NRC(2),HVEKT(ig),F(100),HINT(ig),
     & HBUF(600),JBUF(600),TSFSTA(ig),TSAV(ig),XLFT(NROWXLFT,8,3),
     1 PP(3),TF(14),MNO(3),CP(13),XK(7),PYA(3),POWAB(17,3),
     2 MAKK(3),AIRU(13,3),SFX(13,3), TINT(ig), TVEKT(ig), TBUF(600)
      COMMON/NULL_COM/FACT(nh4),RFACT(nh4),FACTM(nh4),RFACTM(nh4)
     & ,MAA,IFD1,IFD2,KCD,KBCD,NHCD,NHBCD,NNC,IFPL(3),
     & NN1,NN2,NHKTA,NHKTB,NHKTC,NHKTD,KHKTA,KHKTB,KHKTC,KHKTD,NNB
c      COMMON/NULL_COM/FACT(nh4),RFACT(nh4),FACTM(nh4),RFACTM(nh4)
      COMMON /REP/ NEWIND(MXCBF) , MSTOLD(8)
      COMMON /INDX/ PC(512),DSTRT(8,MXCBF),NTAP,LU2,NRSS,NUCZ,ITAG,
     1 MAXLOP,MAXLOT,KMAX,NMAX,KHKT(7),MULT(8),ISYTYP(3),ITYPE(7,28),
     2 AND(8,8),OR(8,8),EOR(8,8),NPARSU(8),NPAR(8),MULNUC(Mxatms),
     3 NHKT(MXTNSH),MUL(MXTNSH),NUCO(MXTNSH),NRCO(MXTNSH),JSTRT(MXTNSH),
     4 NSTRT(MXTNSH),MST(MXTNSH),JRS(MXTNSH)
      COMMON /DAT/  ALPHA(MXTNPR),CONT(MXTNCC),CENT(3,MXTNSH),
     1              CORD(Mxatms,3),CHARGE(Mxatms),FMULT(8),TLA, TLC
      COMMON /HARMON/ HARMT
      CHARACTER*8 TEXT(2)
      CHARACTER*8 TTEXT(2)
      DATA TEXT /'********','ONEHAMIL'/
      DATA TTEXT /'********','KINETINT'/
      DATA LUSCR /13/

#ifdef _TRACE_VMOL
      write(*,*) '@VMOL: Entered file ',__FILE__
#endif

      IPRINT = IFLAGS(1)
C
CCDIR$ BOUNDS
C
      INDECP=0
      IVARV = 0
      REWIND LUSCR
      CP(1)=1.
      PI=2.0*DACOS(0.0D00)
          MINUS1=-1
      NUT=0
      DO 3 I=1,3
    3 POWAB(1,I)=1.
      ISAB=0
      MIA=0
      DO 1 IA=1,KMAX
      NHKTA=NHKT(IA)
      NHKTAM=NHKTA-1
      NANA=(NHKTAM*NHKTA*(NHKTA+1))/6
      KHKTA=KHKT(NHKTA)
      MULA=MUL(IA)
      NUCA=NUCO(IA)
      NUC(1) = NUCA
      JSTA=JSTRT(IA)
      MIB=0
      DO 40 IB=1,IA
       MRAF = NRCO(IA) * NRCO(IB)
       NRC(1) = NRCO(IA)
       NRC(2) = NRCO(IB)
      NHKTB=NHKT(IB)
      NHKTBP=NHKTB+1
      NHKTBQ=NHKTB+2
      INCR=(NHKTA*NHKTBQ*(NHKTA+NHKTBQ))/2
      MAXALL=NHKTA+NHKTB
      NBNB=((NHKTB-1)*NHKTB*NHKTBP)/6
      KHKTB=KHKT(NHKTB)
      MULB=MUL(IB)
      MAB=OR(MULA,MULB)
      KAB=AND(MULA,MULB)
      NUCB=NUCO(IB)
      NUCAB=NUCA*NUCB
      NUC(2) = NUCB
      JJRS(1) = JRS(IA)
      JJRS(2) = JRS(IB)
      JSTB=JSTRT(IB)
      CALL EXCOEF(TSFSTA,NHKTA,NHKTB+2,NUCA,NUCB,JSTA,JSTB,FACTM)
C
C....    SECOND LAST CALL PARAMETER HERE SHOULD BE ENOUGH FOR 8 PRIMITIVE
C....    D'S IN A CONTRACTION.  SEE 2E ROUTINES FOR DEFN
C
CJDW The required parameter is NROWXLFT (we think). It should let us
C    handle larger primitive d sets (like 14) !
c     CALL EXPA (IA,IB,SABB,FACT,XLFT,ALP,BET,S1,CSS1,PV,NPL,CAB,
c    $ 4095,8)
      CALL EXPA (IA,IB,SABB,FACT,XLFT,ALP,BET,S1,CSS1,PV,NPL,CAB,
     $ NROWXLFT,8)
      IF(IA .EQ. IB) CALL BCOMP(SABB,NUCA)
      GO TO 491
      IJAB = 0
      DO 490 N=1,NMAX
      MABN = OR(KAB,MULNUC(N))
      DO 1540 KN=1,MAXLOT
      IF (AND(KN,MABN) .NE. 1) GO TO 1540
      DO 112 KB=1,MAXLOT
      IF (AND(KB,MAB) .NE. 1) GO TO 112
      DO 223 JAB=1,NUCAB
C     THE LOCATION OF THE "T" VALUES MUST BE CHANGED IF
C     THIS SECTION IS TO BE ACTIVATED.
      F(JAB + IJAB) =
     & ((PV(JAB,KB,1) - CORD(N,1)*XAND(1,KN,2))**2 +
     &  (PV(JAB,KB,2) - CORD(N,2)*XAND(2,KN,2))**2 +
     &  (PV(JAB,KB,3) - CORD(N,3)*XAND(3,KN,2))**2 ) * S1(JAB)
  223 CONTINUE
      IJAB = IJAB + NUCAB
  112 CONTINUE
 1540 CONTINUE
  490 CONTINUE
  491 CONTINUE
C     CALL VGAMMA (F,F,MAXALL-2,IJAB)
      MAA=KHKTA*KHKTB
      CALL ZERO (TINT,MRAF*MAA*MULT(KAB))
      CALL ZERO (HINT,MRAF*MAA*MULT(KAB))
      HKAB=FMULT(KAB)
      DO 5 KB=1,MAXLOT
      IF(AND(KB,MAB).NE.1) GO TO 5
      DISTAB=0.
      DO 50 I=1,3
      POWAB(2,I) = CENT(I,IA) - CENT(I,IB)*XAND(I,KB,2)
      DO 1510 J=2,MAXALL
 1510 POWAB(J+1,I) = POWAB(J,I)*POWAB(2,I)
   50 DISTAB = DISTAB+POWAB(3,I)
      INDX=1
      DO 77 I=1,NUCAB
  77  SAAB(I) = CAB(I) * SABB(I,KB)
      JAB = 0
      DO 8 JA=1,NUCA
      JJA=JA+JSTA
      A=ALPHA(JJA)
      DO 8 JB=1,NUCB
      JAB= JAB + 1
      JJB=JB+JSTB
      B=ALPHA(JJB)
      C=A+B
      C2=2.*C
      CINV=1.0/C
      EPSI=0.25*CINV
      EPSJ=0.25*C
      DEL=0.5*CINV
      AT2=1.0
      TF(1)=1.0
      DO 121 N=2,MAXALL
      TF(N)=TF(N-1)*AT2*DEL
  121 AT2=AT2+2.0
      PIC=PI*CINV
      ROPIC=SQRT(PIC)
      TETA = 2./ROPIC
      BCI=B*CINV
      DO 10 I=1,3
      PYA(I) =-BCI*POWAB(2,I)
      PP (I) = CENT(I,IA) + PYA(I)
   10 CONTINUE
      SAABI = PIC*ROPIC*SABB(JAB,KB)
      JNDX=INDX
      DO 12 MM=0,NHKTAM
      DO 12 LL=0,NHKTBP
      LM=(LL+MM)/2
      DO 13 KK=0,LM
      TSAV(JNDX+LL+MM-KK-KK)=TSFSTA(JNDX+KK+KK)*RFACT(KK+1)*EPSJ**KK
   13 CONTINUE
   12 JNDX=JNDX+LL+MM+1
      NUM=0
      DO 100 NA=1,KHKTA
      DO 100 NB=1,KHKTB
      NUM=NUM+1
      DO 56 I=1,3
   56 MNO(I) = NXYZ(NBNB+NB,I) * (NXYZ(NBNB+NB,I)-1)/2
C.....CALCULATE KINETIC ENERGY INTEGRALS ........
      LYDI(1)=0
      LYDI(2)=0
      DO 150 I=1,7
      LYDI(3)=LYDI(2)
      LYDI(2)=LYDI(1)
      LYDI(1)=((I-3*(I/3))/2) * (2-(I/8)) * ((-1)**I)
      SUM=1.0
      DO 140 IJK=1,3
      MM = NXYZ(NANA+NA,IJK)
      LL = NXYZ(NBNB+NB,IJK) + LYDI(IJK)
      IF(LL.LT.0) GO TO 150
      MMM=(LL+MM)/2+1
      IML=INDX+(MM*NHKTBQ*(MM+NHKTBQ)+LL*(MM+MM+LL+1))/2
      ANS=0.
      DO 1101 KK=1,MMM
      ANS=ANS+TSAV(IML+LL+MM-KK-KK+2) *POWAB(3+LL+MM-2*KK,IJK)
 1101 CONTINUE
      SUM=SUM*ANS
      IF(I.GT.1) GO TO 140
      MIJK = NXYZ(NANA+NA,IJK) + NXYZ(NBNB+NB,IJK) +1
      DO 2 KK=1,MIJK
      SFX(KK,IJK) = TSFSTA(IML-1+KK)*POWAB(MIJK-KK+1,IJK)*C**(KK-1)
    2 CONTINUE
  140 CONTINUE
  150 XK(I)=SUM
      XKN=SAABI*(B*((2*NHKTB+1)*XK(1)-2.0*B*(XK(2)+XK(3)+XK(4)))-
     & (MNO(1)*XK(5) + MNO(2)*XK(6) + MNO(3)*XK(7)))
      HVEKT(JAB+(NUM-1)*2*NUCAB) = XKN
      TVEKT(JAB+(NUM-1)*2*NUCAB) = XKN
      TNAI=0.0
C...  CALCULATE NUCLEAR ATTRACTION INTEGRALS
      DO 49 N=1,NMAX
      MU=MULNUC(N)
      MABN=OR(KAB,MU)
      KABN=AND(KAB,MU)
      DO 154 KN=1,MAXLOT
      IF(AND(KN,MABN) .NE. 1) GO TO 154
      T=0.
      DO 301 IJK=1,3
      CP(2) = PP(IJK)-CORD(N,IJK)*XAND(IJK,KN,2)
      CP(3) = 0.5*CP(2)**2
      DO 51 J=5,MAXALL
   51 CP(J-1)=CP(J-2)*CP(2)/FLOAT(J-2)
      T=T+CP(3)
      MAX= NXYZ(NANA+NA,IJK) + NXYZ(NBNB+NB,IJK) +1
      DO 301 I=1,MAX
      AIRU(I,IJK)=0.0
      AI=SFX(I,IJK)
      IRMAX=(I+1)/2
      DO 301 IR=1,IRMAX
      IRUMAX=IRMAX-IR+1
      DO 301 IRU=1,IRUMAX
      IQ=IR+IRU-2
      IP=I-IQ-IQ-1
  301 AIRU(IP+IRU,IJK)=AIRU(IP+IRU,IJK)+AI*CP(IP+1)*EPSI**IQ
     1  *RFACT(IR)/FACTM(IRU)
      INDF=MAX0(MAXALL,6)
      F(INDF)=C2*T
      IONE=1
cold  CALL VGAMMA (F,F,MAXALL-2,IONE)
      CALL VGAMMA2(F,F(INDF),F(INDF+IONE),F(INDF+2*IONE),
     &             MAXALL-2,IONE)
      PETA=TETA*FMULT(KABN)*CHARGE(N)/HKAB
      PN=0.0
      DO 302 I=1,NXYZ(NANA+NA,1) + NXYZ(NBNB+NB,1) +1
      DO 302 J=1,NXYZ(NANA+NA,2) + NXYZ(NBNB+NB,2) +1
      DO 302 K=1,NXYZ(NANA+NA,3) + NXYZ(NBNB+NB,3) +1
  302 PN=PN+AIRU(I,1)*AIRU(J,2)*AIRU(K,3)*F(I+J+K-2)
      TNAI=TNAI-PN*PETA
  154  CONTINUE
   49 CONTINUE
      HVEKT(JAB+(NUM-1)*2*NUCAB) = HVEKT(JAB+(NUM-1)*2*NUCAB)+TNAI*SAABI
  100 CONTINUE
    8 INDX=INDX+INCR
C
      DO 151 I=1,MAA
      CALL TRANS (2,CONT,JJRS,TVEKT(1+(I-1)*2*NUCAB),NUC,NRC,ISTO)
      CALL TRANS (2,CONT,JJRS,HVEKT(1+(I-1)*2*NUCAB),NUC,NRC,ISTO)
      do 151 j=1,mraf
      tvekt((i-1)*mraf+j) = tvekt((i-1)*2*nucab+j)
      hvekt((i-1)*mraf+j) = hvekt((i-1)*2*nucab+j)
  151 CONTINUE
cch---------------------------------------------------------------------
c     ecp: add ecpints to oneints
cch---------------------------------------------------------------------
CSSS      if (iecp.eq.1) then
CSSS        do 152 i=1,maa
CSSS          do 152 j=1,mraf
CSSS            indecp = indecp+1  
cYAU            write(6,*) ecpvect(indecp)
CSSS            hvekt((i-1)*mraf+j) = hvekt((i-1)*mraf+j)+
CSSS     &      ecpvect(indecp)
CSSS 152    continue
CSSS      endif
cch---------------------------------------------------------------------
      JSTO=ISTO
      IF(HARMT) THEN
      CALL HARM1(KHKTA,KHKTB,NHKTA,NHKTB,MAA,MRAF,TVEKT,HLP1)
      CALL HARM1(KHKTA,KHKTB,NHKTA,NHKTB,MAA,MRAF,HVEKT,HLP1)
      ENDIF
C
C...  FORM SYMMETRY INTEGRALS....
      ISYM=0
      DO 18 LA=1,MAXLOP
      NUM=0
      DO 19 NA=1,KHKTA
      IVA=ITYPE(NHKTA,NA)
      NEA=EOR(LA,IVA)
      IF(AND(MULA,NEA).NE.1) GO TO 19
      DO 20 NB=1,KHKTB
      IVA=ITYPE(NHKTB,NB)
      NEB=EOR(LA,IVA)
      IF(AND(MULB,NEB).NE.1) GO TO 20
      ISYM=ISYM+1
      NUN=NUM+NB
      SIGNE=PC(AND(KB,NEB))
      DO 21 KK=1,MRAF
      TINT(KK+(ISYM-1)*MRAF)=TINT(KK+(ISYM-1)*MRAF)+
     + SIGNE*TVEKT(KK+ISTO+(NUN-1)*mraf)
      HINT(KK+(ISYM-1)*MRAF)=HINT(KK+(ISYM-1)*MRAF)+
     + SIGNE*HVEKT(KK+ISTO+(NUN-1)*mraf)
  21  CONTINUE
   20 CONTINUE
   19 NUM=NUM+KHKTB
   18 CONTINUE
    5 ISAB=ISAB+NUCA*NUCB

C
C...  ORDER SYMMETRY INTEGRALS AND ACCOUNT FOR MULTIPLICITY
      ISYM=0
      DO 27 LA=1,MAXLOP
      DO 37 NA=1,KHKTA
      IF(AND(MULA,EOR(LA,ITYPE(NHKTA,NA))) . NE . 1) GO TO 37
      DO 38 NB=1,KHKTB
      IF(AND(MULB,EOR(LA,ITYPE(NHKTB,NB))) . NE . 1) GO TO 38
      ISYM=ISYM+1
      IF(IA.EQ.IB. AND .NA.LT.NB) GO TO 38
      KK=0
      DO 39 K1=1,NRC(1)
      DO 39 K2=1,NRC(2)
      KK=KK+1
      IF(IA.EQ.IB .AND. NA.EQ.NB .AND. K1.LT.K2) GO TO 39
      LBU = DSTRT(LA,NSTRT(IA)+NA) - NRC(1) + K1
      MBU = DSTRT(LA,NSTRT(IB)+NB) - NRC(2) + K2
      LBU = NEWIND(LBU + MSTOLD(LA))
      MBU = NEWIND(MBU + MSTOLD(LA))
      IF(LBU*MBU  .EQ.  0) GO TO 39
      LBU = LBU - MST(LA)
      MBU = MBU - MST(LA)
      NBU = MAX0(LBU,MBU)
      IBUF = NBU*(NBU-3)/2 + LBU + MBU + NPARSU(LA)
      NUT=NUT+1
      TBUF(NUT)=TINT(KK+(ISYM-1)*MRAF)*HKAB
      HBUF(NUT)=HINT(KK+(ISYM-1)*MRAF)*HKAB
      JBUF(NUT)=IBUF
      if(iprint.ge.500)then
       write(6,5000)ibuf,tbuf(nut),hbuf(nut)
5000   format(t5,15x,i5,2f20.12)
      endif
      IF(NUT.EQ.600) THEN
      WRITE (LUSCR) TBUF,HBUF,JBUF,NUT
      IVARV = IVARV + 1
      NUT=0
      ENDIF
   39 CONTINUE
   38 CONTINUE
   37 CONTINUE
   27 CONTINUE
   40 MIB=MIB+KHKTB
    1 MIA=MIA+KHKTA
      IF (NUT.GT.0) THEN
        WRITE (LUSCR) TBUF,HBUF,JBUF,NUT
        IVARV = IVARV + 1
      ENDIF
C
C....    NOW WRITE ONE-ELECTRON INTEGRALS ON NTAP
C
      REWIND LUSCR
      WRITE(NTAP) TEXT(1)
      WRITE(NTAP) TEXT(2)
      DO 600 I = 1, IVARV
      READ(LUSCR) TBUF, HBUF, JBUF, NUT
      WRITE(NTAP) HBUF, JBUF, NUT
600   CONTINUE
#ifdef _DEBUG_VMOL
      print 6655,(hbuf(k), k=1,nut)
 6655 format (10f10.4)
#endif /* _DEBUG_VMOL */
      WRITE(NTAP) HBUF, JBUF, MINUS1
      REWIND LUSCR
      WRITE(NTAP) TTEXT(1)
      WRITE(NTAP) TTEXT(2)
      DO 610 I = 1, IVARV
      READ(LUSCR) TBUF, HBUF, JBUF, NUT
      WRITE(NTAP) TBUF, JBUF, NUT
610   CONTINUE
      WRITE(NTAP) TBUF, JBUF, MINUS1
      REWIND LUSCR

#ifdef _TRACE_VMOL
      write(*,*) '@VMOL: Leaving file ',__FILE__
#endif

      RETURN
      END
