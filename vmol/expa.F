      SUBROUTINE EXPA(IP,IQ,SABB,FAC,XLFT,ALP,BET,S1,CS1,PV,NPL,CCA,
     $ NDIM,MAXLOT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE RFA
      INTEGER DSTRT,AND,OR,EOR

      integer left1,left2,right1,right2,right3,right4

C-----------------------------------------------------------------------
C     Parameters
C-----------------------------------------------------------------------
#include "baslims.par"
#include "mxatms.par"
      parameter (nh4=4*NHT-3, nh2=NHT+NHT+1)
      parameter (mxp2=maxprim*maxprim)
      parameter (nh43=3*nh4,nh4mx=nh4*mxp2)
C-----------------------------------------------------------------------
      LOGICAL NPL(3,8)
      COMMON /FLAG/ IFLAG
      LOGICAL IFLAG
      COMMON /XA/ XAND(3,8,nh4)
      COMMON /INDX/ PC(512),DSTRT(8,MXCBF),NTAP,LU2,NRSS,NUCZ,ITAG,
     1 MAXLOP,MAXXXX,KMAX,NMAX,KHKT(7),MULT(8),ISYTYP(3),ITYPE(7,28),
     2 AND(8,8),OR(8,8),EOR(8,8),NPARSU(8),NPAR(8),MULNUC(Mxatms),
     3 NHKT(MXTNSH),MUL(MXTNSH),NUCO(MXTNSH),NRCO(MXTNSH),JSTRT(MXTNSH),
     4 NSTRT(MXTNSH),MST(MXTNSH),JRS(MXTNSH)
      COMMON /DAT/  ALPHA(MXTNPR),CONT(MXTNCC),CENT(3,MXTNSH),
     1              CORD(Mxatms,3),CHARGE(Mxatms),FMULT(8),TLA, TLC
      DIMENSION SABB(mxp2,8),POW(3,nh4),TSF(10000),XLFT(NDIM,MAXLOT,3),
     1 RFA(nh4),S1(mxp2),GP(mxp2,nh4),SUM(mxp2,3),FAC(nh4),CCA(mxp2),
     2 ALP(mxp2),BET(mxp2),PV(mxp2,8,3),CS1(mxp2),AB(mxp2),ARG(mxp2)

#ifdef _TRACE_VMOL
      write(*,*) '@VMOL: Entered file ',__FILE__
#endif

      do j = 1, nh4
         pow(1,j) = 1.d0
         pow(2,j) = 1.d0
         pow(3,j) = 1.d0
         do i = 1, mxp2
            gp(i,j) = 1.d0
         end do
      end do
C
CCDIR$ BOUNDS
C
      NHKTA=NHKT(IP)
      NHKTB=NHKT(IQ)
      NHP=NHKTA+NHKTB
      NUCA=NUCO(IP)
      NUCB=NUCO(IQ)
      JSTA=JSTRT(IP)
      JSTB=JSTRT(IQ)
      JSA = JRS(IP)
      JSB = JRS(IQ)
      LM=NHP-1
      DO 9 JA=1,NUCA
      DO 9 JB=1,NUCB
      CCA(JB+(JA-1)*NUCB) = CONT(JSA+JA)*CONT(JSB+JB)
      ALP(JB+(JA-1)*NUCB) = ALPHA(JSTA+JA)
      BET(JB+(JA-1)*NUCB) = ALPHA(JSTB+JB)
    9 CONTINUE
      NUCAB=NUCA*NUCB
      DO 1 JPQ=1,NUCAB
      S1(JPQ) = ALP(JPQ) + BET(JPQ)
      CS1(JPQ)=1./S1(JPQ)
      ARG(JPQ) = ALP(JPQ)*BET(JPQ)*CS1(JPQ)
    1 CONTINUE
      DO 8 I=1,3
      DO 8 KQ=1,MAXLOT
      DO 8 JPQ=1,NUCAB
      PV(JPQ,KQ,I)=(CENT(I,IP)*ALP(JPQ) + CENT(I,IQ)*BET(JPQ)*XAND(
     1 I,KQ,2))*CS1(JPQ)
    8 CONTINUE
      DO 12 JAB=1,NUCAB
      AB(JAB)=ALP(JAB)
   12 TSF(JAB)=1.
C
      IM=0
      IMM=IM
      IB=0
      DO 500 L=1,NHKTA
      DO 400 M=1,NHKTB
      LMX=L+M-2
c      IF(LMX-1) 400,300,100

      if(lmx .gt. 1)then
  100   I9 = NUCAB * (M+IB)
        DO 200 K=2,LMX
        I6 = NUCAB*(IM+K-1)
        I7 = I6-I9
CDIR$ IVDEP
*VOCL LOOP,NOVREC
        DO 200 J=1,NUCAB
  200   TSF(J+I6) = (TSF(J+I7+NUCAB)*AB(J) + TSF(J+I7))*CS1(J)
      endif

      if(lmx .ge. 1)then

 300    I0 = NUCAB * (IM+1-M-IB)
        I9 = NUCAB * (IM+LMX)
        I3 = NUCAB * (IM+L-2-IB)
C
C IVDEP DIRECTIVE ON CRAY, NOVECTOR ON FPS 500.
C
CDIR$ IVDEP
*VOCL LOOP,NOVREC
        DO 310 J=1,NUCAB
        TSF(J + NUCAB*IM) = TSF(J + I0)*AB(J)*CS1(J)
  310   TSF(J+I9) = TSF(J+I3)*CS1(J)
      endif

  400 IM=IM+LMX+1
      IF (L.EQ.1) THEN
      DO 410 JAB=1,NUCAB
  410 AB(JAB) = -BET(JAB)
      ENDIF
  500 IB=KHKT(NHKTB)+NHKTB*(L-1)
C
      DO 610 L=1,NHKTA
      DO 610 M=1,NHKTB
      LMX=L+M-1
      DO 600 ML=1,LMX
      IMM=IMM+1
      call sscal(nucab,fac(ml),tsf(1+nucab*(imm-1)),1)
  600 CONTINUE
  610 CONTINUE
C
C FPS COMPILER CHOKES ON THIS!
C      DO 600 JAB=1,NUCAB
C  600 TSF(JAB + NUCAB*(IMM-1))=TSF(JAB + NUCAB*(IMM-1))*FAC(ML)
C
      DO 11 I=3,LM,2
      DO 11 J=1,NUCAB
   11 GP(J,I) = 0.5*S1(J)*GP(J,I-2)/DFLOAT(I-1)

      DO 900 KB=1,MAXLOT

      DIST=0.
      DO 830 K=1,3
      POW(K,2)=CENT(K,IP)-CENT(K,IQ)*XAND(K,KB,2)
      POW(K,3)=POW(K,2)*POW(K,2)
      NPL(K,KB) = POW(K,3) .LT. 1.E-12
C     IF(NPL(K,KB)) GO TO 444
C     IF(NPL(K,KB)) GO TO 830
      IF(NHP .GE. 5)THEN
        DO 820 IND=5,NHP
        POW(K,IND-1)=POW(K,IND-2)*POW(K,2)
  820   CONTINUE
      ENDIF
      DIST=DIST+POW(K,3)
  830 CONTINUE
C
      DO 840 J=1,NUCAB
      ARGJ = MIN(ARG(J)*DIST,80.D0)
      SABB(J,KB)= EXP(-ARGJ)
  840 CONTINUE

      DO 890 I=1,3
Cold  NUM=0
      DO 880 NC=1,NHKTA
      DO 880 ND=1,NHKTB
      I34=NC+ND-1
C
C     New way to calculate NUM+I34. This does not depend on any
C     recursion. Probably doesn't help, but you never know, and at
C     least I know what is going on. NUMALT replaces NUM+I34 in the
C     calculation of I4.
C
      NUMALT = (NHKTB*(NC-1)*NC)/2 - (NC-1)*NHKTB +
     &         ((NC-1)*NHKTB*(NHKTB+1))/2 +
     &         (ND*(ND+1))/2 +
     &          ND*(NC-1)
C
      INCR = 1
C     IF(NPL(I,KB)) INCR = 2
C
      DO 870 IND=1,I34,INCR
      CALL ZERO(SUM(1,I),NUCAB)
      IR2M = 1
C     IF(NPL(I,KB)) IR2M = IND
Cold  I4 = NUCAB*(NUM+I34-IND)
      I4 = NUCAB*(NUMALT-IND)
C
C     Replace these two loops by subroutine expa_inner to see this
C     routine will work with default optimization level at OSC.
C
      CALL EXPA_INNER(SUM,GP,TSF,POW,MXP2,NH4,I,I4,IR2M,IND,NUCAB)
ctest      DO 860 IR2=IR2M,IND,2
ctest      I5 = I4+NUCAB*(IR2-1)
ctest      DO 850 J=1,NUCAB
ctest      SUM(J,I) = SUM(J,I) + GP(J,IR2)*TSF(J+I5)*POW(I,IND-IR2+1)
ctest  850 CONTINUE
ctest  860 CONTINUE
C
c      DO 40 J=1,NUCAB
c   40 XLFT(J+I4,KB,I) = SUM(J,I)*RFA(I34-IND+1)
C
      scale=rfa(i34-ind+1)
      call scopy(nucab,sum(1,i),1,xlft(1+i4,kb,i),1)
      call sscal(nucab,scale,xlft(1+i4,kb,i), 1)

  870 CONTINUE
Cold  NUM = NUM + I34
  880 CONTINUE
  890 CONTINUE
  900 CONTINUE
C
      IF(IP.EQ.IQ) CALL ACOMP(SABB,NUCA)
      RETURN
      ENTRY RFINIT
      RFA(1)=1.
      do 151 i=1,nh4-1
151   rfa(i+1)=0.25*rfa(i)/dfloat(i)

#ifdef _TRACE_VMOL
      write(*,*) '@VMOL: Leaving file ',__FILE__
#endif

      return
      END
