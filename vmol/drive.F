      SUBROUTINE DRIVE (I2,A2,I1,A1,IA,IC,ID,NPL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DSTRT,AND,OR,EOR
C-----------------------------------------------------------------------
C     Parameters
C-----------------------------------------------------------------------
#include "baslims.par"
#include "mxatms.par"
      parameter (nh4=4*nht-3, nh2=nht+nht+1)
      parameter (khtt=nht*(nht+1)/2, kh4=khtt**4)
      parameter (mxp2=maxprim*maxprim)
C-----------------------------------------------------------------------
      LOGICAL NPL(3,8,2),IFLAG,IFPL,IFD1,IFD2,IFD1X,IFD2X,HARMT
      DIMENSION I2(20),A2(1),I1(50),A1(1)
      COMMON /FLAG/ IFLAG
      COMMON /TAB/ NNNGAM
      COMMON /REP/ NEWIND(MXCBF) , MSTOLD(8)
      COMMON /ITY/ ITYA(28),ITYB(28),ITYC(28),ITYD(28),
     1 KPQRS,NSTA,NSTB,NSTC,NSTD,MUT,NUMSM(4),NBLOCK(4),
     1 JA(8),JB(8),JC(8),JD(8),BUF(2400),IBUF(2400),NUTTY(4),N2TAPE(4),
     X LEOR8(64),IFD3,MULA,MULB,MULC,MULD
      COMMON /MMMM/ M2,M3,LINA,MUM,MIM,NMAS,MX,MRAF,NMSS,LIA,JJUNK(2)
      COMMON /XA/ XAND(3,8,nh4)
      COMMON /DAT/  ALPHA(MXTNPR),CONT(MXTNCC),CENT(3,MXTNSH),
     1              CORD(Mxatms,3),CHARGE(Mxatms),FMULT(8),TLA, TLC
      COMMON /ONE/ 
     1 ALP(mxp2),BET(mxp2),S1(mxp2),CSS1(mxp2),PV(mxp2,8,3),CAB(mxp2),
     2 GAM(mxp2),DEL(mxp2),S2(mxp2),CSS2(mxp2),QV(mxp2,8,3),CCD(mxp2)
      COMMON/NULL_COM/FACT(nh4),RFACT(nh4),FACTM(nh4),RFACTM(nh4)
     & ,MAA,IFD1,IFD2,KCD,KBCD,NHCD,NHBCD,NNC,IFPL(3),
     & NN1,NN2,NHKTA,NHKTB,NHKTC,NHKTD,KHKTA,KHKTB,KHKTC,KHKTD,NNB
      COMMON /CON/ RAX,RAY,RAZ,RPRX,RPRY,RPRZ,DISTPQ,DISTRS,TOLR,TP52,
     1 SABB(mxp2,8),SCDD(mxp2,8)
      COMMON /CONI/ NUCR,NUCS,NUCA,NUCQ,
     & NRCR,NRCS,NRCA,NRCQ,JSTR,JSTS,JSTA,JSTQ,
     & IFD1X,IFD2X,NUCAQ,NUCRS,NAQRS,NCQRS,KQ,KR,KS,IVA
      COMMON /INDX/ PC(512),DSTRT(8,MXCBF),NTAP,LU2,NRSS,NUCZ,ITAG,
     1 MAXLOP,MAXLOT,KMAX,NMAX,KHKT(7),MULT(8),ISYTYP(3),ITYPE(7,28),
     2 IVAND(8,8),IVOR(8,8),IVEOR(8,8),NPARSU(8),NPAR(8),MULNUC(Mxatms),
     3 NHKT(MXTNSH),MUL(MXTNSH),NUCO(MXTNSH),NRCO(MXTNSH),JSTRT(MXTNSH),
     4 NSTRT(MXTNSH),MST(MXTNSH),JRS(MXTNSH)
c      COMMON /INDX/ PC(512),DSTRT(8,500),NTAP,LU2,NRSS,NUCZ,ITAG,
c     1 MAXLOP,MAXLOT,KMAX,NMAX,KHKT(7),MULT(8),ISYTYP(3),ITYPE(7,28)
      COMMON /HARMON/ HARMT
c      DIMENSION HLP(kh4),HLP2(kh4)
c       DIMENSION HLP(kh4),HLP2(1)
C...  These dimensions are not rigorous!!
C     space requirement is MRAF*MAA.
C

#ifdef _TRACE_VMOL
      write(*,*) '@VMOL: Entered file ',__FILE__
#endif

C...  HIGHER QUANTUM NUMBER CASE
      NUCRS=NUCR*NUCS
      NAQRS = NUCAQ*NUCRS
      I2(1)=1
      I1(4)=I1(3)
      I1(5)=I1(4)
      I1(6)=I1(5)+NAQRS
      I1(7)=I1(6)+NAQRS
      I2(2)=I2(1)+NAQRS
      I2(3)=I2(2)
      I2(4)=I2(3)+NAQRS
      I2(5)=I2(4)+NAQRS
      I2(6)=I2(5)+NAQRS
      KPQ=AND(MULA-1,MULB-1)
      MPQ=OR(MULA-1,MULB-1)
      KPQR=AND(KPQ,MULC-1)
      MPQR=OR(KPQ,MULC-1)
      MPQRS=OR(KPQR,MULD-1)
      MX = AND(KPQR,MULD-1)
C...  BIT"J" OF MX IS 1  IF ALL CENTERS ARE ON THE SYMMETRY ELEMENT "J"
      KPQRS = MX+1
      NMSS = MULT(MULA)*MULT(MULB)*MULT(MULC)*MULT(MULD)
      NMAS = NMSS/MULT(KPQRS)
      NNB = NHKTA*NHBCD*(NNNGAM+1)/2
      NND = NUCAQ*NHKTA*NHKTB*(NN1+1)/2
      NNE = NUCRS*NHCD*(NN2+1)/2
      IF(NMSS .EQ. 1) THEN
C...  NO SYMMETRY
      IFD1X = IFD1
      IFD2X = IFD2
      CALL LOOP1 (A1(I1(3)),A2(I2(2)),I2,A2,I1,A1)
      CALL EXPA(IC,ID,SCDD,FACTM,A1(I1(2)),GAM,DEL,S2,CSS2,QV,
     $ NPL(1,1,2),CCD,NNE,MAXLOT)
      CALL SET2 (NRSS,ITAG,I2,A2,I1,A1)
      DO 71 J=1,3
   71 IFPL(J) = NPL(J,1,1) .AND. NPL(J,1,2) .AND.
     1 ABS(CENT(J,IA) - CENT(J,IC)) .LT. 1.E-07
      IF(ITAG .EQ. 1) THEN
C     ONLY ONE PASS
      CALL SAVR (A1(I1(5)),NAQRS)
      I2(4)=I2(3)+MAA*MRAF
      I2(5)=I2(4)+NAQRS
      I2(6)=I2(5)+NAQRS
      CALL CONLOR (A1(I1(7)),SABB,SCDD,NPL,A2(I2(3)),A1(I1(1)),
     $ A1(I1(2)),A2(I2(4)),A1(I1(5)),A1(I1(6)),A2(I2(5)),NND,
     $ NNE,MAXLOT,I2,A2,I1,A1)
C
      ELSE
C     SEVERAL PASSES
      CALL ZERO(A2(I2(3)),MRAF*MAA)
      I2(4)=I2(3)+MRAF*MAA
      DO 1013 NUCZ=1,NUCR*NUCS,NRSS
      NUCRS = NRSS
      NAQRS = NUCAQ*NUCRS
      CALL SAVVR(A1(I1(5)),S2(NUCZ),CSS2(NUCZ),NAQRS)
      I2(5)=I2(4)+NAQRS
      I2(6)=I2(5)+NAQRS
      I1(6)=I1(5)+NAQRS
      I1(7)=I1(6)+NAQRS
      CALL CONLOX (A1(I1(7)),SABB,SCDD(NUCZ,1),NPL,A2(I2(3)),
     $ A1(I1(1)),A1(I1(2)+NUCZ-1),A2(I2(4)),A1(I1(5)),
     $ A1(I1(6)),A2(I2(5)),QV(NUCZ,1,1),NND,NNE,MAXLOT,I2,A2,I1,A1)
 1013 CONTINUE
      ENDIF
      IF(.NOT.IFD1X.AND.NRCA.GT.1)CALL AAEXP(A2(I2(3)),NRCA,NRCR*NRCS)
      IF(.NOT.IFD2X.AND.NRCR.GT.1)CALL BBEXP(A2(I2(3)),NRCR,NRCA*NRCQ)
C
      IF(HARMT)
c     $ CALL HARMONY(MRAF,MAA*MRAF,A2(I2(3)),HLP,HLP2)
     $ CALL HARMONY(MRAF,MAA*MRAF,A2(I2(3)),A2(I2(4)),ZJUNK,MAA)
C
      CALL LOOPN3 (A1(I1(3)),A2(I2(3)),I2,A2,I1,A1)
C
      ELSE
C...  SYMMETRY
      CALL LOOP1 (A1(I1(3)),A2(I2(2)),I2,A2,I1,A1)
      CALL ZERO (A1(I1(4)),LINA*MRAF)
      CALL EXPA(IC,ID,SCDD,FACTM,A1(I1(2)),GAM,DEL,S2,CSS2,QV,
     $ NPL(1,1,2),CCD,NNE,MAXLOT)
      CALL SET2(NRSS,ITAG,I2,A2,I1,A1)
      NUCRS=NRSS
      NAQRS=NUCAQ*NUCRS
      I2(4)=I2(3)+MAA*MRAF
      I2(5)=I2(4)+NAQRS
      I2(6)=I2(5)+NAQRS
      I1(6)=I1(5)+NAQRS
      I1(7)=I1(6)+NAQRS
C...   LOOP OVER SYMMETRY TRANSFORMATIONS
      DO 84 KQ=1,MAXLOP
      IFD1X = IFD1 .OR. (KQ.GT.1)
      IF(AND(KQ-1,MPQ  ) .NE.0) GO TO 84
      DO 83 KR=1,MAXLOP
      IF(AND(KR-1,MPQR ) .NE.0) GO TO 83
      DO 82 KS=1,MAXLOP
      IF(AND(KS-1,MPQRS) .NE.0) GO TO 82
      IVA = NEQV(KR-1,KS-1) + 1
      IFD2X=IFD2.OR.KR.NE.KS
      DO 81 J=1,3
   81 IFPL(J) = NPL(J,KQ,1) .AND. NPL(J,IVA,2) .AND.
     1 ABS(CENT(J,IA) - CENT(J,IC)*XAND(J,KR,2)) .LT. 1.E-07
      CALL ZERO(A2(I2(3)),MRAF*MAA)
      DO 1012 NUCZ=1,NUCR*NUCS,NRSS
      CALL SAVVR(A1(I1(5)),S2(NUCZ),CSS2(NUCZ),NAQRS)
      CALL CONLOX(A1(I1(7)),SABB(1,KQ),SCDD(NUCZ,IVA),NPL,A2(I2(3)),
     $ A1(I1(1)),A1(I1(2)+NUCZ-1),A2(I2(4)),A1(I1(5)),A1(I1(6)),
     $ A2(I2(5)),QV(NUCZ,1,1),NND,NNE,MAXLOT,I2,A2,I1,A1)
 1012 CONTINUE
      IF(.NOT.IFD1X.AND.NRCA.GT.1)CALL AAEXP(A2(I2(3)),NRCA,NRCR*NRCS)
      IF(.NOT.IFD2X.AND.NRCR.GT.1)CALL BBEXP(A2(I2(3)),NRCR,NRCA*NRCQ)
C
      IF(HARMT)
c     $ CALL HARMONY(MRAF,MAA*MRAF,A2(I2(3)),HLP,HLP2)
     $ CALL HARMONY(MRAF,MAA*MRAF,A2(I2(3)),A2(I2(4)),ZJUNK,MAA)
C
      K=((KQ-1)*MAXLOP + KR-1)*MAXLOP + KS-1
      CALL LMUL (A1(I1(4)),A2(I2(3)),A1(I1(3)),MRAF,NMAS,LIA,A2(I2(2))
     $ ,K,A1(I1(7)),TLC)
   82 CONTINUE
   83 CONTINUE
   84 CONTINUE
cOLD  CALL LOOP3 (A1(I1(3)),A1(I1(4)),I2,A2,I1,A1)
      CALL LOOP3 (A1(I1(3)),A1(I1(4)))
      ENDIF

#ifdef _TRACE_VMOL
      write(*,*) '@VMOL: Leaving file ',__FILE__
#endif

      RETURN
      END
