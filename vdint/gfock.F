

      SUBROUTINE GFOCK(FHFA,FHFB,NF,ISOR,ISOS,ISOT,
     &                 NORBA,NORBB,NORBC,NORBD,
     &                 NOABCD,ICGTOA,ICGTOB,ICGTOC,ICGTOD,
     &                 DIAGAB,DIAGCD,SHABAB,SHAEQB,SHCEQD,
     &                 TCONAB,TCONCD,
     &                 KHKTA,KHKTB,KHKTC,KHKTD,
     &                 NHKTA,NHKTB,NHKTC,NHKTD,
     &                 THRESH,SYMFAC,IPRINT,
     &                 WORK1,LWORK1,
     &                 DHFA,DHFB,NBASIS,
     &                 NCOORD,NIR)
C
C     CONSTRUCTS THE DERIVATIVES OF THE FOCK MATRICES WITH
C     RESPECT TO THE COMPONENTS OF AN EXTERNAL MAGNETIC FIELD
C
CEND
C 
C   SEP/91 JG UNIVERSITY OF KARLSRUHE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
      LOGICAL TWOCEN, THRCEN, FOUCEN, DERONE, DERTWO
      LOGICAL NODER, DIAGAB, DIAGCD, SHABAB, SHAEQB, SHCEQD,
     &        BIGVEC, DCMPAB, DCMPCD, TCONAB, TCONCD, NOPMAT, NODV
      LOGICAL SCF,NONHF
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
#include "mxcent.par"
      DIMENSION FOCK1(3),FOCK2(3),
     &          FOCK3(3),FOCK4(3),FOCK5(3),
     &          FOCK6(3),FOCK3B(3),FOCK4B(3),FOCK5B(3),
     &          FOCK6B(3)
      DIMENSION ICGTOA(MXAOVC), ICGTOB(MXAOVC),
     &          ICGTOC(MXAOVC), ICGTOD(MXAOVC)
      DIMENSION WORK1(LWORK1)
      DIMENSION DHFA(NBASIS*(NBASIS+1)/2,0:NIR),
     &          DHFB(NBASIS*(NBASIS+1)/2,0:NIR)
      DIMENSION FHFA(NBASIS*(NBASIS+1)/2,0:NIR,NCOORD),
     &          FHFB(NBASIS*(NBASIS+1)/2,0:NIR,NCOORD)
C
      COMMON /INTADR/ IWKAO, IWKSO, IWKHHS, IWK1HH, IWK1HC, IWKLST
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     &                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     &                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     &                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     &                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     &                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     &                IPTXYZ(3,0:7)
      COMMON/METHOD/IUHF,SCF,NONHF
      COMMON/EXPCOM/SIGNXYZ(12),
     &              NCENT1, NCENT2, NCENT3, NCENT4,ISON(4),
     &              DERONE, DERTWO, TWOCEN, THRCEN, FOUCEN,
     &              NINTYP, NCCINT
      COMMON/PERT/NTPERT,NPERT(8),IPERT(8),IXPERT,IYPERT,IZPERT,
     &            IXYZPERT(3)
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
C  FACTORS :   A.EQ.B        FACTOR OF HALF (EXCEPT DIAGONAL)
C  FACTORS :   C.EQ.D        FACTOR OF HALF  (EXCEPT DIAGONAL)
C  FACTORS :    AB EQ. CD    FACTOR OF HALF 
C
      SFAC=SYMFAC
      IF(SHABAB) SFAC=SFAC*HALF
      IF(SHAEQB) SFAC=SFAC*HALF
      IF(SHCEQD) SFAC=SFAC*HALF
C
      THRSH=MAX(THRESH,ZERADD)
C
C FIRST SET UP LOOP FOR RHF
C 
      IF(IUHF.EQ.0) THEN
C
C  Loop over shell components 
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
C  Loop over contracted functions 
C
          DO 500 IORBA=1,NORBA
           INDXA=ICGTOA(IORBA)+IORBA+(ICOMPA-1)*NORBA
           NORBBB=NORBB
           IF (TCONAB) NORBBB=IORBA
           DO 600 IORBB=1,NORBBB
            INDXB=ICGTOB(IORBB)+IORBB+(ICOMPB-1)*NORBB
            DHFAB=DHFA(ITRI(INDXA,INDXB),ISOR)*HALF
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
              DHFCD=DHFA(ITRI(INDXC,INDXD),ISOS)*HALF
              FCD=ONE
              IF(TCONCD.AND.IORBC.NE.IORBD) FCD=TWO
               NODER=.FALSE.
               INT1=IAOADR
               INT2=INT1+3*NCCINT  
               FACT=FCABAB*FAB*FCD
C
C LOOP OVER X, Y, AND Z COMPONENT OF MAGNETIC FIELD
C
               DO 900 I = 1, 3
C
C CALCULATE THE TWO COULOMB CONTRIBUTIONS
C
                FOCK1(I)=DHFAB*(WORK1(INT1)-WORK1(INT2))*FACT
                FOCK2(I)=DHFCD*(WORK1(INT1)+WORK1(INT2))*FACT
c                write(*,*) 'integral 1 ',work1(int1)
c                write(*,*) 'integral 2 ',work1(int2)
c                write(*,*) i,fock1(i),INDXC,INDXD,DHFAB
c                write(*,*) i,fock2(i),INDXA,INDXB,DHFCD
C
C CALCULATE THE FOUR EXCHANGE CONTRIBUTIONS
C
                FOCK3(I)=-FOURTH*DHFAC*WORK1(INT2)*FACT
                FOCK4(I)=-FOURTH*DHFBD*WORK1(INT2)*FACT
                FOCK5(I)=-FOURTH*DHFAD*WORK1(INT1)*FACT
                FOCK6(I)=-FOURTH*DHFBC*WORK1(INT1)*FACT
C
c                write(*,*) i,fock3(i),INDXB,INDXD,DHFAC
c                write(*,*) i,fock4(i),INDXA,INDXC,DHFBD
c                write(*,*) i,fock5(i),INDXC,INDXB,DHFAD
c                write(*,*) i,fock6(i),INDXA,INDXD,DHFBC
                INT1=INT1+NCCINT
                INT2=INT2+NCCINT
900            CONTINUE
               IREP=0
C
C DISTRIBUTE NOW THE FOCK1 AND FOCK2 ELEMENTS TO THE REAL FOCK MATRIX 
C
              DO 904 I=1,3
C
               IXYZ=IXYZPERT(I)-1
               XB=PT(IBTAND(IXYZ,ISOR))
               XC=PT(IBTAND(IXYZ,ISOT))
               XD=PT(IBTAND(IXYZ,ISOTS))
C
C WHAT IS THE CORRECT SIGN ?
C
               IF(INDXC.LE.INDXD) THEN
               FHFA(ITRI(INDXC,INDXD),ISOS,I)=
     &               FHFA(ITRI(INDXC,INDXD),ISOS,I)+FOCK1(I)
     &               *XC
c       write(*,*) 'fock',fhfa(itri(INDXC,INDXD),0,i),xc,INDXC,INDXD
               ENDIF
               IF(INDXC.GE.INDXD) THEN
               FHFA(ITRI(INDXC,INDXD),ISOS,I)=
     &               FHFA(ITRI(INDXC,INDXD),ISOS,I)-FOCK1(I)
     &               *XD
c       write(*,*) 'fock',fhfa(itri(INDXC,INDXD),0,i),xd,INDXC,INDXD
               ENDIF
C
               IF(INDXA.LE.INDXB) THEN
               FHFA(ITRI(INDXA,INDXB),ISOR,I)=
     &               FHFA(ITRI(INDXA,INDXB),ISOR,I)+FOCK2(I)
c       write(*,*) 'fock',fhfa(itri(INDXA,INDXB),0,i),one,INDXA,INDXB
               ENDIF
               IF(INDXA.GE.INDXB) THEN
               FHFA(ITRI(INDXA,INDXB),ISOR,I)=
     &               FHFA(ITRI(INDXA,INDXB),ISOR,I)-FOCK2(I)
     &               *XB
c       write(*,*) 'fock',fhfa(itri(INDXA,INDXB),0,i),xb,INDXA,INDXB
               ENDIF
C
C EXCHANGE CONTRIBUTIONS
C
               IF(INDXB.LE.INDXD) THEN                     
               FHFA(ITRI(INDXB,INDXD),ISORTS,I)=
     &               FHFA(ITRI(INDXB,INDXD),ISORTS,I)-FOCK3(I)
     &               *XB
c       write(*,*) 'fock',fhfa(itri(INDXB,INDXD),0,i),INDXB,INDXD
               ENDIF
               IF(INDXB.GE.INDXD) THEN
               FHFA(ITRI(INDXB,INDXD),ISORTS,I)=
     &               FHFA(ITRI(INDXB,INDXD),ISORTS,I)+FOCK3(I)
     &               *XD
c       write(*,*) 'fock',fhfa(itri(INDXB,INDXD),0,i),INDXB,INDXD
               ENDIF
               IF(INDXA.LE.INDXC) THEN
               FHFA(ITRI(INDXA,INDXC),ISOT,I)=
     &               FHFA(ITRI(INDXA,INDXC),ISOT,I)+FOCK4(I)
c       write(*,*) 'fock',fhfa(itri(INDXA,INDXC),0,i),INDXA,INDXC
               ENDIF
               IF(INDXA.GE.INDXC) THEN
               FHFA(ITRI(INDXA,INDXC),ISOT,I)=
     &               FHFA(ITRI(INDXA,INDXC),ISOT,I)-FOCK4(I)
     &               *XC
c       write(*,*) 'fock',fhfa(itri(INDXA,INDXC),0,i),INDXA,INDXC
               ENDIF
               IF(INDXB.LE.INDXC) THEN
               FHFA(ITRI(INDXB,INDXC),ISORT,I)=
     &               FHFA(ITRI(INDXB,INDXC),ISORT,I)-FOCK5(I)
     &               *XB
c       write(*,*) 'fock',fhfa(itri(INDXB,INDXC),0,i),INDXB,INDXC
               ENDIF
               IF(INDXB.GE.INDXC) THEN
               FHFA(ITRI(INDXB,INDXC),ISORT,I)=
     &               FHFA(ITRI(INDXB,INDXC),ISORT,I)+FOCK5(I)
     &               *XC
c       write(*,*) 'fock',fhfa(itri(INDXB,INDXC),0,i),INDXB,INDXC
               ENDIF
               IF(INDXA.LE.INDXD) THEN
               FHFA(ITRI(INDXA,INDXD),ISOTS,I)=
     &               FHFA(ITRI(INDXA,INDXD),ISOTS,I)+FOCK6(I)
c       write(*,*) 'fock',fhfa(itri(INDXA,INDXD),0,i),INDXA,INDXD
               ENDIF
               IF(INDXA.GE.INDXD) THEN 
                FHFA(ITRI(INDXA,INDXD),ISOTS,I)=
     &               FHFA(ITRI(INDXA,INDXD),ISOTS,I)-FOCK6(I)
     &               *XD
c       write(*,*) 'fock',fhfa(itri(INDXA,INDXD),0,i),INDXA,INDXD
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
C  LOOP OVER SHELL COMPONENTS
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
C  LOOP OVER CONTRACTED FUNCTIONS 
C
          DO 2500 IORBA=1,NORBA
           INDXA=ICGTOA(IORBA)+IORBA+(ICOMPA-1)*NORBA
           NORBBB=NORBB
           IF (TCONAB) NORBBB=IORBA
           DO 2600 IORBB=1,NORBBB
            INDXB=ICGTOB(IORBB)+IORBB+(ICOMPB-1)*NORBB
            DHFAB=(DHFA(ITRI(INDXA,INDXB),ISOR)+
     &             DHFB(ITRI(INDXA,INDXB),ISOR))*HALF
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
              DHFCD=(DHFA(ITRI(INDXC,INDXD),ISOS)+
     &              DHFB(ITRI(INDXC,INDXD),ISOS))*HALF
              FCD=ONE
              IF(TCONCD.AND.IORBC.NE.IORBD) FCD=TWO
               NODER = .FALSE.
               INT1= IAOADR
               INT2=INT1+3*NCCINT
               FACT=FCABAB*FAB*FCD
               DO 2900 I = 1, 3 
C
C CALCULATE THE TWO COLOUMB-CONTRIBUTIONS
C
                FOCK1(I)=DHFAB*(WORK1(INT1)+WORK1(INT2))*FACT
                FOCK2(I)=DHFCD*(WORK1(INT1)-WORK1(INT2))*FACT
C
C CALCULATE THE EIGHT EXCHANGE CONTRIBUTIONS
C
                FOCK3(I)=-HALF*DHFAAC*WORK1(INT2)*FACT
                FOCK3B(I)=-HALF*DHFBAC*WORK1(INT2)*FACT
                FOCK4(I)=-HALF*DHFABD*WORK1(INT2)*FACT
                FOCK4B(I)=-HALF*DHFBBD*WORK1(INT2)*FACT
                FOCK5(I)=-HALF*DHFAAD*WORK1(INT1)*FACT
                FOCK5B(I)=-HALF*DHFBAD*WORK1(INT1)*FACT
                FOCK6(I)=-HALF*DHFABC*WORK1(INT1)*FACT
                FOCK6B(I)=-HALF*DHFBBC*WORK1(INT1)*FACT
C
                INT1=INT1+NCCINT
                INT2=INT2+NCCINT
2900           CONTINUE
C
          IREP=0
C
C DISTRIBUTE NOW THE FOCK1 AND FOCK2 ELEMENTS TO THE REAL FOCK MATRIX 
C
              DO 2904 I=1,3
C
               IXYZ=IXYZPERT(I)-1
               XB=PT(IBTAND(IXYZ,ISOR))
               XC=PT(IBTAND(IXYZ,ISOT))
               XD=PT(IBTAND(IXYZ,ISOTS))
C
               IF(INDXC.LE.INDXD)THEN
                FHFA(ITRI(INDXC,INDXD),ISOS,I)=
     &               FHFA(ITRI(INDXC,INDXD),ISOS,I)+FOCK1(I)*
     &               XC
                FHFB(ITRI(INDXC,INDXD),ISOS,I)=
     &               FHFB(ITRI(INDXC,INDXD),ISOS,I)+FOCK1(I)*
     &               XC
               ENDIF
               IF(INDXC.GE.INDXD)THEN
                FHFA(ITRI(INDXC,INDXD),ISOS,I)=
     &               FHFA(ITRI(INDXC,INDXD),ISOS,I)-FOCK1(I)*
     &               XD
                FHFB(ITRI(INDXC,INDXD),ISOS,I)=
     &               FHFB(ITRI(INDXC,INDXD),ISOS,I)-FOCK1(I)*
     &               XD
               ENDIF
               IF(INDXA.LE.INDXB)THEN
                FHFA(ITRI(INDXA,INDXB),ISOR,I)=
     &               FHFA(ITRI(INDXA,INDXB),ISOR,I)+FOCK2(I)
                FHFB(ITRI(INDXA,INDXB),ISOR,I)=
     &               FHFB(ITRI(INDXA,INDXB),ISOR,I)+FOCK2(I)
               ENDIF
               IF(INDXA.GE.INDXB)THEN
                FHFA(ITRI(INDXA,INDXB),ISOR,I)=
     &               FHFA(ITRI(INDXA,INDXB),ISOR,I)-FOCK2(I)*
     &               XB
                FHFB(ITRI(INDXA,INDXB),ISOR,I)=
     &               FHFB(ITRI(INDXA,INDXB),ISOR,I)-FOCK2(I)*
     &               XB
               ENDIF
               IF(INDXB.LE.INDXD)THEN
                FHFA(ITRI(INDXB,INDXD),ISORTS,I)=
     &               FHFA(ITRI(INDXB,INDXD),ISORTS,I)-FOCK3(I)*
     &               XB
                FHFB(ITRI(INDXB,INDXD),ISORTS,I)=
     &               FHFB(ITRI(INDXB,INDXD),ISORTS,I)-FOCK3B(I)*
     &               XB
               ENDIF
               IF(INDXB.GE.INDXD)THEN
                FHFA(ITRI(INDXB,INDXD),ISORTS,I)=
     &               FHFA(ITRI(INDXB,INDXD),ISORTS,I)+FOCK3(I)*
     &               XD
                FHFB(ITRI(INDXB,INDXD),ISORTS,I)=
     &               FHFB(ITRI(INDXB,INDXD),ISORTS,I)+FOCK3B(I)*
     &               XD
               ENDIF
               IF(INDXA.LE.INDXC)THEN
                FHFA(ITRI(INDXA,INDXC),ISOT,I)=
     &               FHFA(ITRI(INDXA,INDXC),ISOT,I)+FOCK4(I)
                FHFB(ITRI(INDXA,INDXC),ISOT,I)=
     &               FHFB(ITRI(INDXA,INDXC),ISOT,I)+FOCK4B(I)
               ENDIF 
               IF(INDXA.GE.INDXC)THEN
                FHFA(ITRI(INDXA,INDXC),ISOT,I)=
     &               FHFA(ITRI(INDXA,INDXC),ISOT,I)-FOCK4(I)*
     &               XC
                FHFB(ITRI(INDXA,INDXC),ISOT,I)=
     &               FHFB(ITRI(INDXA,INDXC),ISOT,I)-FOCK4B(I)*
     &               XC
               ENDIF
               IF(INDXB.LE.INDXC)THEN
                FHFA(ITRI(INDXB,INDXC),ISORT,I)=
     &               FHFA(ITRI(INDXB,INDXC),ISORT,I)-FOCK5(I)*
     &               XB
                FHFB(ITRI(INDXB,INDXC),ISORT,I)=
     &               FHFB(ITRI(INDXB,INDXC),ISORT,I)-FOCK5B(I)*
     &               XB
               ENDIF
               IF(INDXB.GE.INDXC)THEN
                FHFA(ITRI(INDXB,INDXC),ISORT,I)=
     &               FHFA(ITRI(INDXB,INDXC),ISORT,I)+FOCK5(I)*
     &               XC
                FHFB(ITRI(INDXB,INDXC),ISORT,I)=
     &               FHFB(ITRI(INDXB,INDXC),ISORT,I)+FOCK5B(I)*
     &               XC
               ENDIF
               IF(INDXA.LE.INDXD)THEN
                FHFA(ITRI(INDXA,INDXD),ISOTS,I)=
     &               FHFA(ITRI(INDXA,INDXD),ISOTS,I)+FOCK6(I)
                FHFB(ITRI(INDXA,INDXD),ISOTS,I)=
     &               FHFB(ITRI(INDXA,INDXD),ISOTS,I)+FOCK6B(I)
               ENDIF
               IF(INDXA.GE.INDXD)THEN
                FHFA(ITRI(INDXA,INDXD),ISOTS,I)=
     &               FHFA(ITRI(INDXA,INDXD),ISOTS,I)-FOCK6(I)*
     &               XD
                FHFB(ITRI(INDXA,INDXD),ISOTS,I)=
     &               FHFB(ITRI(INDXA,INDXD),ISOTS,I)-FOCK6B(I)*
     &               XD
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
