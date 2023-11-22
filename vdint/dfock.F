      SUBROUTINE DFOCK(FHFA,FHFB,NF,ISOR,ISOS,ISOT,
     *                 NORBA,NORBB,NORBC,NORBD,
     *                 NOABCD,ICGTOA,ICGTOB,ICGTOC,ICGTOD,
     *                 DIAGAB,DIAGCD,SHABAB,SHAEQB,SHCEQD,
     *                 TCONAB,TCONCD,
     *                 KHKTA,KHKTB,KHKTC,KHKTD,
     *                 NHKTA,NHKTB,NHKTC,NHKTD,
     *                 THRESH,SYMFAC,IPRINT,
     *                 WORK1,LWORK1,
     *                 DHFA,DHFB,NBASIS,
     *                 NCOORD,NIR)
C
C     CONSTRUCTS THE DERIVATIVES OF THE FOCK MATRICES WITH
C     RESPECT TO THE NUCLEAR COORDINATES
C
CEND
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
      LOGICAL TWOCEN, THRCEN, FOUCEN, DERONE, DERTWO
      LOGICAL NODER, DIAGAB, DIAGCD, SHABAB, SHAEQB, SHCEQD,
     *        BIGVEC, DCMPAB, DCMPCD, TCONAB, TCONCD, NOPMAT, NODV
      LOGICAL SCF,NONHF,hfdftgrad
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
#include "mxcent.par"
      DIMENSION IPERT(12),FOCK1(12),FOCK2(12),
     *          FOCK3(12),FOCK4(12),FOCK5(12),
     *          FOCK6(12),FOCK3B(12),FOCK4B(12),FOCK5B(12),
     *          FOCK6B(12)
      DIMENSION ICGTOA(MXAOVC), ICGTOB(MXAOVC),
     *          ICGTOC(MXAOVC), ICGTOD(MXAOVC)
C
      COMMON /INTADR/ IWKAO, IWKSO, IWKHHS, IWK1HH, IWK1HC, IWKLST
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
      DIMENSION WORK1(LWORK1)
C
      DIMENSION DHFA(NBASIS*(NBASIS+1)/2,0:NIR),
     *          DHFB(NBASIS*(NBASIS+1)/2,0:NIR)
      DIMENSION FHFA(NBASIS*(NBASIS+1)/2,0:NIR,NCOORD),
     *          FHFB(NBASIS*(NBASIS+1)/2,0:NIR,NCOORD)
C
C  NEW COMMON BLOCK IN THE ACES II ENVIRONMENT FOR KEEPINT TRACK ABOUT
C  THE METHOD 
C
      COMMON/METHOD/IUHF,SCF,NONHF
C
      COMMON /EXPCOM/  SIGNXYZ(12),
C     *                SIGN1X, SIGN1Y, SIGN1Z, SIGN2X, SIGN2Y, SIGN2Z,
C     *                SIGN3X, SIGN3Y, SIGN3Z, SIGN4X, SIGN4Y, SIGN4Z,
     *                NCENT1, NCENT2, NCENT3, NCENT4,
     *                ISON(4),
c     *                ISO1,   ISO2,   ISO3,   ISO4,
     *                DERONE, DERTWO, TWOCEN, THRCEN, FOUCEN,
     *                NINTYP, NCCINT
C
      COMMON /FLAGS2/ IFLAGS2(500)
C
      DATA FOURTH,HALF,ONE,TWO,FOUR,ZERADD,AZERO /0.25D0,0.5D0,1.0D0,
     &     2.0D0,4.0D0,1.D-16,0.D0/
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTXOR(I,J) = XOR(I,J)
C
      ITRI(I,J) = MAX(I,J)*(MAX(I,J) - 1)/2 + MIN(I,J)
      INDEX(I,J,K,L) = NORBA*(NORBB*(NORBC*(L-1)+K-1)+J-1)+I
C
      ISOTS=IBTXOR(ISOT,ISOS)
      ISORT=IBTXOR(ISOR,ISOT)
      ISORTS=IBTXOR(ISOR,ISOTS)
C
C  FACTORS :   A.NE.B.NE.C.NE.D ONE
C  FACTORS :   A.EQ.B    FACTOR OF HALF (EXCEPT DIAGONAL)
C  FACTORS :   C.EQ.D    FACTOR OF HALF  (EXCEPT DIAGONAL)
C  FACTORS :    AB EQ. CD FACPRS OF HALF 
C
      SFAC=SYMFAC
      IF(SHABAB) SFAC=SFAC*HALF
      IF(SHAEQB) SFAC=SFAC*HALF
      IF(SHCEQD) SFAC=SFAC*HALF
C
c      write(*,*) shabab,shaeqb,shceqd,sfac,tconab,tconcd 
C
C  NUMBER OF FIRST DERIVATIVE INTEGRALS
C
      IF(TWOCEN) NUMDER=3
      IF(THRCEN) NUMDER=6
      IF(FOUCEN) NUMDER=9
C
      THRSH=MAX(THRESH,ZERADD)
C
C For HF-DFT gradients. P. Verma and A. Perera, 08/2008
C
       hfdftgrad=.false.
       call getrec(-1,'JOBARC','HFDFTGRA',1,ihfdftgrad)
       if (ihfdftgrad .eq. 1) Then
           IFLAGS2(153) = 0
           hfdftgrad    = .True.
       endif
C 
C FIRST SET UP LOOP FOR RHF
C
      IF(IUHF.EQ.0) THEN
C
C     ***** Loop over shell components *****
C
      IAOFF=IWKAO+1
      DO 100 ICOMPA=1,KHKTA
       KHKTBB=KHKTB
       IF(DIAGAB) KHKTBB=ICOMPA
       DO 200 ICOMPB = 1,KHKTBB
        DCMPAB=DIAGAB.AND.ICOMPA.EQ.ICOMPB
        FACAB=PT(IBTAND(ITYPE(NHKTB,ICOMPB),ISOR))
        IF(DIAGAB.AND.ICOMPA.NE.ICOMPB) FACAB=TWO*FACAB
        DO 300 ICOMPC=1,KHKTC
         PTC=PT(IBTAND(ITYPE(NHKTC,ICOMPC),ISOT))
         KHKTDD=KHKTD
         IF(DIAGCD) KHKTDD=ICOMPC
         DO 400 ICOMPD=1,KHKTDD
          DCMPCD=DIAGCD.AND.ICOMPC.EQ.ICOMPD
          FACCD=PTC*PT(IBTAND(ITYPE(NHKTD,ICOMPD),ISOTS))
          IF(DIAGCD.AND.ICOMPC.NE.ICOMPD) FACCD=TWO*FACCD
          FCABAB=FACAB*FACCD*SFAC
          IAOADR=IAOFF
C
C  Loop over contracted functions *****
C
          DO 500 IORBA=1,NORBA
           INDXA=ICGTOA(IORBA)+IORBA+(ICOMPA-1)*NORBA
           NORBBB=NORBB
           IF (TCONAB) NORBBB=IORBA
           DO 600 IORBB=1,NORBBB
            INDXB=ICGTOB(IORBB)+IORBB+(ICOMPB-1)*NORBB
            DHFAB=DHFA(ITRI(INDXA,INDXB),ISOR)
            FAB=ONE
            IF(TCONAB.AND.IORBA.NE.IORBB) FAB=TWO
            DO 700 IORBC=1,NORBC
             INDXC=ICGTOC(IORBC)+IORBC+(ICOMPC-1)*NORBC
             DHFAC=DHFA(ITRI(INDXA,INDXC),ISOT)
             DHFBC=DHFA(ITRI(INDXB,INDXC),ISORT)
             NORBDD=NORBD
             IF(TCONCD) NORBDD=IORBC
             DO 800 IORBD=1,NORBDD
              INDXD=ICGTOD(IORBD)+IORBD+(ICOMPD-1)*NORBD
              DHFAD=DHFA(ITRI(INDXA,INDXD),ISOTS)
              DHFBD=DHFA(ITRI(INDXB,INDXD),ISORTS)
              DHFCD=DHFA(ITRI(INDXC,INDXD),ISOS)
              FCD=ONE
              IF(TCONCD.AND.IORBC.NE.IORBD) FCD=TWO
               NODER = .FALSE.
               INT   = IAOADR
#ifdef _DEBUG
               if (int+NCCINT*NUMDER.gt.lwork1) then
                  print *, '@DFOCK: Out-of-bounds read.'
                  call errex
               end if
#endif
               FACT=FCABAB*FAB*FCD
c               write(*,*)indxa,indxb,indxc,indxd,fact
               DO 900 I = 1, NUMDER
C
C CALCULATE THE TWO COULOMB CONTRIBUTIONS
                FOCK1(I)=DHFAB*WORK1(INT)*FACT
                FOCK2(I)=DHFCD*WORK1(INT)*FACT
C CALCULATE THE FOUR EXCHANGE CONTRIBUTIONS
                FOCK3(I)=-FOURTH*DHFAC*WORK1(INT)*FACT
                FOCK4(I)=-FOURTH*DHFBD*WORK1(INT)*FACT
                FOCK5(I)=-FOURTH*DHFAD*WORK1(INT)*FACT
                FOCK6(I)=-FOURTH*DHFBC*WORK1(INT)*FACT
C
                INT = INT + NCCINT
900            CONTINUE
C
C CALCULATE THE ``MISSING'' INTEGRAL DERIVATIVE PIECE
C
               DO 899 IREP=0,MAXLOP
               IF(TWOCEN) THEN
                IPERT(1)=IPTCNT(3*NCENT1-2,IREP)
                IPERT(2)=IPTCNT(3*NCENT1-1,IREP)
                IPERT(3)=IPTCNT(3*NCENT1,IREP)
                IPERT(4)=IPTCNT(3*NCENT2-2,IREP)
                IPERT(5)=IPTCNT(3*NCENT2-1,IREP)
                IPERT(6)=IPTCNT(3*NCENT2,IREP)
                DO 901 I=1,3
                 FOCK1(I+3)=-FOCK1(I)
                 FOCK2(I+3)=-FOCK2(I)
                 FOCK3(I+3)=-FOCK3(I)
                 FOCK4(I+3)=-FOCK4(I)
                 FOCK5(I+3)=-FOCK5(I)
                 FOCK6(I+3)=-FOCK6(I)
901             CONTINUE
               ELSE IF(THRCEN) THEN
                IPERT(1)=IPTCNT(3*NCENT1-2,IREP)
                IPERT(2)=IPTCNT(3*NCENT1-1,IREP)
                IPERT(3)=IPTCNT(3*NCENT1,IREP)
                IPERT(4)=IPTCNT(3*NCENT2-2,IREP)
                IPERT(5)=IPTCNT(3*NCENT2-1,IREP)
                IPERT(6)=IPTCNT(3*NCENT2,IREP)
                IPERT(7)=IPTCNT(3*NCENT3-2,IREP)
                IPERT(8)=IPTCNT(3*NCENT3-1,IREP)
                IPERT(9)=IPTCNT(3*NCENT3,IREP)
                DO 902 I=1,3
                 FOCK1(I+6)=-(FOCK1(I)+FOCK1(I+3))
                 FOCK2(I+6)=-(FOCK2(I)+FOCK2(I+3))
                 FOCK3(I+6)=-(FOCK3(I)+FOCK3(I+3))
                 FOCK4(I+6)=-(FOCK4(I)+FOCK4(I+3))
                 FOCK5(I+6)=-(FOCK5(I)+FOCK5(I+3))
                 FOCK6(I+6)=-(FOCK6(I)+FOCK6(I+3))
902             CONTINUE
               ELSE IF(FOUCEN) THEN
                IPERT(1)=IPTCNT(3*NCENT1-2,IREP)
                IPERT(2)=IPTCNT(3*NCENT1-1,IREP)
                IPERT(3)=IPTCNT(3*NCENT1,IREP)
                IPERT(4)=IPTCNT(3*NCENT2-2,IREP)
                IPERT(5)=IPTCNT(3*NCENT2-1,IREP)
                IPERT(6)=IPTCNT(3*NCENT2,IREP)
                IPERT(7)=IPTCNT(3*NCENT3-2,IREP)
                IPERT(8)=IPTCNT(3*NCENT3-1,IREP)
                IPERT(9)=IPTCNT(3*NCENT3,IREP)
                IPERT(10)=IPTCNT(3*NCENT4-2,IREP)
                IPERT(11)=IPTCNT(3*NCENT4-1,IREP)
                IPERT(12)=IPTCNT(3*NCENT4,IREP)
                DO 903 I=1,3
                 FOCK1(I+9)=-(FOCK1(I)+FOCK1(I+3)+FOCK1(I+6))
                 FOCK2(I+9)=-(FOCK2(I)+FOCK2(I+3)+FOCK2(I+6))
                 FOCK3(I+9)=-(FOCK3(I)+FOCK3(I+3)+FOCK3(I+6))
                 FOCK4(I+9)=-(FOCK4(I)+FOCK4(I+3)+FOCK4(I+6))
                 FOCK5(I+9)=-(FOCK5(I)+FOCK5(I+3)+FOCK5(I+6))
                 FOCK6(I+9)=-(FOCK6(I)+FOCK6(I+3)+FOCK6(I+6))
903             CONTINUE
               ENDIF 
               NUMTOT=NUMDER+3
               XA=1.0
               XB=PT(IBTAND(IREP,ISOR))
               XC=PT(IBTAND(IREP,ISOT))
               XD=PT(IBTAND(IREP,ISOTS))
c          call zero(fock1,12)
c          call zero(fock2,12)
c          call zero(fock3,12)
c          call zero(fock4,12)
c          call zero(fock5,12)
c          call zero(fock6,12)
C
C DISTRIBUTE NOW THE FOCK1 AND FOCK2 ELEMENTS TO THE REAL FOCK MATRIX 
C
              DO 904 I=1,NUMTOT
               J=(I-1)/3+1
               JXYZ=I-(J-1)*3
               IF(IPERT(I).NE.0) THEN
                XFAC=PT(IBTAND(IREP,ISON(J)))*SIGNXYZ(I)
C
C WHAT IS THE CORRECT SIGN ?
C
                if(indxc.le.indxd) then
                 FHFA(ITRI(INDXC,INDXD),ISOS,IPERT(I))=
     &           FHFA(ITRI(INDXC,INDXD),ISOS,IPERT(I))+FOCK1(I)*XFAC*XC
                endif
                if(indxc.ge.indxd)then
                 FHFA(ITRI(INDXC,INDXD),ISOS,IPERT(I))=
     &           FHFA(ITRI(INDXC,INDXD),ISOS,IPERT(I))+FOCK1(I)*XFAC*XD
                endif
                if(indxa.le.indxb) then
                 FHFA(ITRI(INDXA,INDXB),ISOR,IPERT(I))=
     &           FHFA(ITRI(INDXA,INDXB),ISOR,IPERT(I))+FOCK2(I)*XFAC
                endif
                if(indxa.ge.indxb) then
                 FHFA(ITRI(INDXA,INDXB),ISOR,IPERT(I))=
     &           FHFA(ITRI(INDXA,INDXB),ISOR,IPERT(I))+FOCK2(I)*XFAC*XB
                endif
C
C EXCHANGE BITS
C
                IF (IFLAGS2(153).EQ.0) THEN
                 if(indxb.le.indxd) then                     
                  FHFA(ITRI(INDXB,INDXD),ISORTS,IPERT(I))=
     &            FHFA(ITRI(INDXB,INDXD),ISORTS,IPERT(I))
     &            +FOCK3(I)*XFAC*XB
                 endif
                 if(indxb.ge.indxd) then
                  FHFA(ITRI(INDXB,INDXD),ISORTS,IPERT(I))=
     &            FHFA(ITRI(INDXB,INDXD),ISORTS,IPERT(I))
     &            +FOCK3(I)*XFAC*XD
                 endif
                 if(indxa.le.indxc) then
                  FHFA(ITRI(INDXA,INDXC),ISOT,IPERT(I))=
     &            FHFA(ITRI(INDXA,INDXC),ISOT,IPERT(I))
     &            +FOCK4(I)*XFAC
                 endif
                 if(indxa.ge.indxc) then
                  FHFA(ITRI(INDXA,INDXC),ISOT,IPERT(I))=
     &            FHFA(ITRI(INDXA,INDXC),ISOT,IPERT(I))
     &            +FOCK4(I)*XFAC*XC
                 endif
                 if(indxb.le.indxc) then
                  FHFA(ITRI(INDXB,INDXC),ISORT,IPERT(I))=
     &            FHFA(ITRI(INDXB,INDXC),ISORT,IPERT(I))
     &            +FOCK5(I)*XFAC*XB
                 endif
                 if(indxb.ge.indxc) then
                  FHFA(ITRI(INDXB,INDXC),ISORT,IPERT(I))=
     &            FHFA(ITRI(INDXB,INDXC),ISORT,IPERT(I))
     &            +FOCK5(I)*XFAC*XC
                 endif
                 if(indxa.le.indxd) then
                  FHFA(ITRI(INDXA,INDXD),ISOTS,IPERT(I))=
     &            FHFA(ITRI(INDXA,INDXD),ISOTS,IPERT(I))
     &            +FOCK6(I)*XFAC
                 endif
                 if(indxa.ge.indxd) then
                  FHFA(ITRI(INDXA,INDXD),ISOTS,IPERT(I))=
     &            FHFA(ITRI(INDXA,INDXD),ISOTS,IPERT(I))
     &            +FOCK6(I)*XFAC*XD
                 endif
                ENDIF
               ENDIF
904           CONTINUE
899           CONTINUE
              IAOADR=IAOADR+1
800          CONTINUE
700         CONTINUE
600        CONTINUE
500       CONTINUE
          IAOFF=IAOFF+NOABCD
400      CONTINUE
300     CONTINUE
200    CONTINUE
100   CONTINUE
C
      ELSE
C
C  LOOP FOR UHF
C
C     ***** Loop over shell components *****
C
      IAOFF=IWKAO+1
      DO 2100 ICOMPA=1,KHKTA
       KHKTBB=KHKTB
       IF(DIAGAB) KHKTBB=ICOMPA
       DO 2200 ICOMPB = 1,KHKTBB
        DCMPAB=DIAGAB.AND.ICOMPA.EQ.ICOMPB
        FACAB=PT(IBTAND(ITYPE(NHKTB,ICOMPB),ISOR))
        IF(DIAGAB.AND.ICOMPA.NE.ICOMPB) FACAB=TWO*FACAB
        DO 2300 ICOMPC=1,KHKTC
         PTC=PT(IBTAND(ITYPE(NHKTC,ICOMPC),ISOT))
         KHKTDD=KHKTD
         IF(DIAGCD) KHKTDD=ICOMPC
         DO 2400 ICOMPD=1,KHKTDD
          DCMPCD=DIAGCD.AND.ICOMPC.EQ.ICOMPD
          FACCD=PTC*PT(IBTAND(ITYPE(NHKTD,ICOMPD),ISOTS))
          IF(DIAGCD.AND.ICOMPC.NE.ICOMPD) FACCD=TWO*FACCD
          FCABAB=FACAB*FACCD*SFAC
          IAOADR=IAOFF
C
C  Loop over contracted functions *****
C
          DO 2500 IORBA=1,NORBA
           INDXA=ICGTOA(IORBA)+IORBA+(ICOMPA-1)*NORBA
           NORBBB=NORBB
           IF (TCONAB) NORBBB=IORBA
           DO 2600 IORBB=1,NORBBB
            INDXB=ICGTOB(IORBB)+IORBB+(ICOMPB-1)*NORBB
            DHFAB=DHFA(ITRI(INDXA,INDXB),ISOR)+
     &             DHFB(ITRI(INDXA,INDXB),ISOR)
            FAB=ONE
            IF(TCONAB.AND.IORBA.NE.IORBB) FAB=TWO
            DO 2700 IORBC=1,NORBC
             INDXC=ICGTOC(IORBC)+IORBC+(ICOMPC-1)*NORBC
             DHFAAC=DHFA(ITRI(INDXA,INDXC),ISOT)
             DHFBAC=DHFB(ITRI(INDXA,INDXC),ISOT)
             DHFABC=DHFA(ITRI(INDXB,INDXC),ISORT)
             DHFBBC=DHFB(ITRI(INDXB,INDXC),ISORT)
             NORBDD=NORBD
             IF(TCONCD) NORBDD=IORBC
             DO 2800 IORBD=1,NORBDD
              INDXD=ICGTOD(IORBD)+IORBD+(ICOMPD-1)*NORBD
              DHFAAD=DHFA(ITRI(INDXA,INDXD),ISOTS)
              DHFBAD=DHFB(ITRI(INDXA,INDXD),ISOTS)
              DHFABD=DHFA(ITRI(INDXB,INDXD),ISORTS)
              DHFBBD=DHFB(ITRI(INDXB,INDXD),ISORTS)
              DHFCD=DHFA(ITRI(INDXC,INDXD),ISOS)+
     &              DHFB(ITRI(INDXC,INDXD),ISOS)
              FCD=ONE
              IF(TCONCD.AND.IORBC.NE.IORBD) FCD=TWO
               NODER = .FALSE.
               INT   = IAOADR
#ifdef _DEBUG
               if (int+NCCINT*NUMDER.gt.lwork1) then
                  print *, '@DFOCK: Out-of-bounds read.'
                  call errex
               end if
#endif
               FACT=FCABAB*FAB*FCD
c               write(*,*)indxa,indxb,indxc,indxd,fact
               DO 2900 I = 1, NUMDER
C
C CALCULATE THE TWO COLOUMB-CONTRIBUTIONS
                FOCK1(I)=DHFAB*WORK1(INT)*FACT
                FOCK2(I)=DHFCD*WORK1(INT)*FACT
C CALCULATE THE EIGHT EXCHANGE CONTRIBUTIONS
                FOCK3(I)=-HALF*DHFAAC*WORK1(INT)*FACT
                FOCK3B(I)=-HALF*DHFBAC*WORK1(INT)*FACT
                FOCK4(I)=-HALF*DHFABD*WORK1(INT)*FACT
                FOCK4B(I)=-HALF*DHFBBD*WORK1(INT)*FACT
                FOCK5(I)=-HALF*DHFAAD*WORK1(INT)*FACT
                FOCK5B(I)=-HALF*DHFBAD*WORK1(INT)*FACT
                FOCK6(I)=-HALF*DHFABC*WORK1(INT)*FACT
                FOCK6B(I)=-HALF*DHFBBC*WORK1(INT)*FACT
C
                INT = INT + NCCINT
2900           CONTINUE
C
C CALCULATE THE ``MISSING'' INTEGRAL DERIVATIVE PIECE
C
             DO 2899 IREP=0,MAXLOP
              IF(TWOCEN) THEN
               IPERT(1)=IPTCNT(3*NCENT1-2,IREP)
               IPERT(2)=IPTCNT(3*NCENT1-1,IREP)
               IPERT(3)=IPTCNT(3*NCENT1,IREP)
               IPERT(4)=IPTCNT(3*NCENT2-2,IREP)
               IPERT(5)=IPTCNT(3*NCENT2-1,IREP)
               IPERT(6)=IPTCNT(3*NCENT2,IREP)
               DO 2901 I=1,3
                FOCK1(I+3)=-FOCK1(I)
                FOCK2(I+3)=-FOCK2(I)
                FOCK3(I+3)=-FOCK3(I)
                FOCK4(I+3)=-FOCK4(I)
                FOCK5(I+3)=-FOCK5(I)
                FOCK6(I+3)=-FOCK6(I)
                FOCK3B(I+3)=-FOCK3B(I)
                FOCK4B(I+3)=-FOCK4B(I)
                FOCK5B(I+3)=-FOCK5B(I)
                FOCK6B(I+3)=-FOCK6B(I)
2901           CONTINUE
              ELSE IF(THRCEN) THEN
               IPERT(1)=IPTCNT(3*NCENT1-2,IREP)
               IPERT(2)=IPTCNT(3*NCENT1-1,IREP)
               IPERT(3)=IPTCNT(3*NCENT1,IREP)
               IPERT(4)=IPTCNT(3*NCENT2-2,IREP)
               IPERT(5)=IPTCNT(3*NCENT2-1,IREP)
               IPERT(6)=IPTCNT(3*NCENT2,IREP)
               IPERT(7)=IPTCNT(3*NCENT3-2,IREP)
               IPERT(8)=IPTCNT(3*NCENT3-1,IREP)
               IPERT(9)=IPTCNT(3*NCENT3,IREP)
               DO 2902 I=1,3
                FOCK1(I+6)=-(FOCK1(I)+FOCK1(I+3))
                FOCK2(I+6)=-(FOCK2(I)+FOCK2(I+3))
                FOCK3(I+6)=-(FOCK3(I)+FOCK3(I+3))
                FOCK4(I+6)=-(FOCK4(I)+FOCK4(I+3))
                FOCK5(I+6)=-(FOCK5(I)+FOCK5(I+3))
                FOCK6(I+6)=-(FOCK6(I)+FOCK6(I+3))
                FOCK3B(I+6)=-(FOCK3B(I)+FOCK3B(I+3))
                FOCK4B(I+6)=-(FOCK4B(I)+FOCK4B(I+3))
                FOCK5B(I+6)=-(FOCK5B(I)+FOCK5B(I+3))
                FOCK6B(I+6)=-(FOCK6B(I)+FOCK6B(I+3))
2902           CONTINUE
              ELSE IF(FOUCEN) THEN
               IPERT(1)=IPTCNT(3*NCENT1-2,IREP)
               IPERT(2)=IPTCNT(3*NCENT1-1,IREP)
               IPERT(3)=IPTCNT(3*NCENT1,IREP)
               IPERT(4)=IPTCNT(3*NCENT2-2,IREP)
               IPERT(5)=IPTCNT(3*NCENT2-1,IREP)
               IPERT(6)=IPTCNT(3*NCENT2,IREP)
               IPERT(7)=IPTCNT(3*NCENT3-2,IREP)
               IPERT(8)=IPTCNT(3*NCENT3-1,IREP)
               IPERT(9)=IPTCNT(3*NCENT3,IREP)
               IPERT(10)=IPTCNT(3*NCENT4-2,IREP)
               IPERT(11)=IPTCNT(3*NCENT4-1,IREP)
               IPERT(12)=IPTCNT(3*NCENT4,IREP)
               DO 2903 I=1,3
                FOCK1(I+9)=-(FOCK1(I)+FOCK1(I+3)+FOCK1(I+6))
                FOCK2(I+9)=-(FOCK2(I)+FOCK2(I+3)+FOCK2(I+6))
                FOCK3(I+9)=-(FOCK3(I)+FOCK3(I+3)+FOCK3(I+6))
                FOCK4(I+9)=-(FOCK4(I)+FOCK4(I+3)+FOCK4(I+6))
                FOCK5(I+9)=-(FOCK5(I)+FOCK5(I+3)+FOCK5(I+6))
                FOCK6(I+9)=-(FOCK6(I)+FOCK6(I+3)+FOCK6(I+6))
                FOCK3B(I+9)=-(FOCK3B(I)+FOCK3B(I+3)+FOCK3B(I+6))
                FOCK4B(I+9)=-(FOCK4B(I)+FOCK4B(I+3)+FOCK4B(I+6))
                FOCK5B(I+9)=-(FOCK5B(I)+FOCK5B(I+3)+FOCK5B(I+6))
                FOCK6B(I+9)=-(FOCK6B(I)+FOCK6B(I+3)+FOCK6B(I+6))
2903           CONTINUE
              ENDIF 
               NUMTOT=NUMDER+3
               XA=1.0
               XB=PT(IBTAND(IREP,ISOR))
               XC=PT(IBTAND(IREP,ISOT))
               XD=PT(IBTAND(IREP,ISOTS))
C
C DISTRIBUTE NOW THE FOCK1 AND FOCK2 ELEMENTS TO THE REAL FOCK MATRIX 
C
              DO 2904 I=1,NUMTOT
               J=(I-1)/3+1 
               IF(IPERT(I).NE.0) THEN
C
                XFAC=PT(IBTAND(IREP,ISON(J)))*SIGNXYZ(I)
C
                IF(INDXC.LE.INDXD)THEN
                 FHFA(ITRI(INDXC,INDXD),ISOS,IPERT(I))=
     &           FHFA(ITRI(INDXC,INDXD),ISOS,IPERT(I))+FOCK1(I)*XFAC*XC
                 FHFB(ITRI(INDXC,INDXD),ISOS,IPERT(I))=
     &           FHFB(ITRI(INDXC,INDXD),ISOS,IPERT(I))+FOCK1(I)*XFAC*XC
                ENDIF
                IF(INDXC.GE.INDXD)THEN
                 FHFA(ITRI(INDXC,INDXD),ISOS,IPERT(I))=
     &           FHFA(ITRI(INDXC,INDXD),ISOS,IPERT(I))+FOCK1(I)*XFAC*XD
                 FHFB(ITRI(INDXC,INDXD),ISOS,IPERT(I))=
     &           FHFB(ITRI(INDXC,INDXD),ISOS,IPERT(I))+FOCK1(I)*XFAC*XD
                ENDIF
                IF(INDXA.LE.INDXB)THEN
                 FHFA(ITRI(INDXA,INDXB),ISOR,IPERT(I))=
     &           FHFA(ITRI(INDXA,INDXB),ISOR,IPERT(I))+FOCK2(I)*XFAC
                 FHFB(ITRI(INDXA,INDXB),ISOR,IPERT(I))=
     &           FHFB(ITRI(INDXA,INDXB),ISOR,IPERT(I))+FOCK2(I)*XFAC
                ENDIF
                IF(INDXA.GE.INDXB)THEN
                 FHFA(ITRI(INDXA,INDXB),ISOR,IPERT(I))=
     &           FHFA(ITRI(INDXA,INDXB),ISOR,IPERT(I))+FOCK2(I)*XFAC*XB
                 FHFB(ITRI(INDXA,INDXB),ISOR,IPERT(I))=
     &           FHFB(ITRI(INDXA,INDXB),ISOR,IPERT(I))+FOCK2(I)*XFAC*XB
                ENDIF
                IF (IFLAGS2(153).EQ.0) THEN
                 IF(INDXB.LE.INDXD)THEN
                  FHFA(ITRI(INDXB,INDXD),ISORTS,IPERT(I))=
     &            FHFA(ITRI(INDXB,INDXD),ISORTS,IPERT(I))
     &            +FOCK3(I)*XFAC*XB
                  FHFB(ITRI(INDXB,INDXD),ISORTS,IPERT(I))=
     &            FHFB(ITRI(INDXB,INDXD),ISORTS,IPERT(I))
     &            +FOCK3B(I)*XFAC*XB
                 ENDIF
                 IF(INDXB.GE.INDXD)THEN
                  FHFA(ITRI(INDXB,INDXD),ISORTS,IPERT(I))=
     &            FHFA(ITRI(INDXB,INDXD),ISORTS,IPERT(I))
     &            +FOCK3(I)*XFAC*XD
                  FHFB(ITRI(INDXB,INDXD),ISORTS,IPERT(I))=
     &            FHFB(ITRI(INDXB,INDXD),ISORTS,IPERT(I))
     &            +FOCK3B(I)*XFAC*XD
                 ENDIF
                 IF(INDXA.LE.INDXC)THEN
                  FHFA(ITRI(INDXA,INDXC),ISOT,IPERT(I))=
     &            FHFA(ITRI(INDXA,INDXC),ISOT,IPERT(I))
     &            +FOCK4(I)*XFAC
                  FHFB(ITRI(INDXA,INDXC),ISOT,IPERT(I))=
     &            FHFB(ITRI(INDXA,INDXC),ISOT,IPERT(I))
     &            +FOCK4B(I)*XFAC
                 ENDIF
                 IF(INDXA.GE.INDXC)THEN
                  FHFA(ITRI(INDXA,INDXC),ISOT,IPERT(I))=
     &            FHFA(ITRI(INDXA,INDXC),ISOT,IPERT(I))
     &            +FOCK4(I)*XFAC*XC
                  FHFB(ITRI(INDXA,INDXC),ISOT,IPERT(I))=
     &            FHFB(ITRI(INDXA,INDXC),ISOT,IPERT(I))
     &            +FOCK4B(I)*XFAC*XC
                 ENDIF
                 IF(INDXB.LE.INDXC)THEN
                  FHFA(ITRI(INDXB,INDXC),ISORT,IPERT(I))=
     &            FHFA(ITRI(INDXB,INDXC),ISORT,IPERT(I))
     &            +FOCK5(I)*XFAC*XB
                  FHFB(ITRI(INDXB,INDXC),ISORT,IPERT(I))=
     &            FHFB(ITRI(INDXB,INDXC),ISORT,IPERT(I))
     &            +FOCK5B(I)*XFAC*XB
                 ENDIF
                 IF(INDXB.GE.INDXC)THEN
                  FHFA(ITRI(INDXB,INDXC),ISORT,IPERT(I))=
     &            FHFA(ITRI(INDXB,INDXC),ISORT,IPERT(I))
     &            +FOCK5(I)*XFAC*XC
                  FHFB(ITRI(INDXB,INDXC),ISORT,IPERT(I))=
     &            FHFB(ITRI(INDXB,INDXC),ISORT,IPERT(I))
     &            +FOCK5B(I)*XFAC*XC
                 ENDIF
                 IF(INDXA.LE.INDXD)THEN
                  FHFA(ITRI(INDXA,INDXD),ISOTS,IPERT(I))=
     &            FHFA(ITRI(INDXA,INDXD),ISOTS,IPERT(I))
     &            +FOCK6(I)*XFAC
                  FHFB(ITRI(INDXA,INDXD),ISOTS,IPERT(I))=
     &            FHFB(ITRI(INDXA,INDXD),ISOTS,IPERT(I))
     &            +FOCK6B(I)*XFAC
                 ENDIF
                 IF(INDXA.GE.INDXD)THEN
                  FHFA(ITRI(INDXA,INDXD),ISOTS,IPERT(I))=
     &            FHFA(ITRI(INDXA,INDXD),ISOTS,IPERT(I))
     &            +FOCK6(I)*XFAC*XD
                  FHFB(ITRI(INDXA,INDXD),ISOTS,IPERT(I))=
     &            FHFB(ITRI(INDXA,INDXD),ISOTS,IPERT(I))
     &            +FOCK6B(I)*XFAC*XD
                 ENDIF
                ENDIF
               ENDIF
2904          CONTINUE
2899         CONTINUE
             IAOADR=IAOADR+1
2800         CONTINUE
2700        CONTINUE
2600       CONTINUE
2500      CONTINUE
          IAOFF=IAOFF+NOABCD
2400     CONTINUE
2300    CONTINUE
2200   CONTINUE
2100  CONTINUE
C
      ENDIF
C
      RETURN
      END
