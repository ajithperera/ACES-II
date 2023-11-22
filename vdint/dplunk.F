      SUBROUTINE DPLUNK(SO,ISOR,ISOS,ISOT,
     &                 NORBA,NORBB,NORBC,NORBD,
     &                 NOABCD,ICGTOA,ICGTOB,ICGTOC,ICGTOD,
     &                 DIAGAB,DIAGCD,TCONAB,TCONCD,
     &                 KHKTA,KHKTB,KHKTC,KHKTD,
     &                 NHKTA,NHKTB,NHKTC,NHKTD,
     &                 THRESH,SYMFAC,IPRINT,WORK1,LWORK1,
     &                 NBASIS,NCOORD,NIR,
     &                 BUCK,IBUCK,IBKSIZ,NINBCK,ICHAIN,NREC,
     &                 IMAP)

C
C
C WRITE INTEGRAL DERIVATIVES TO SORT FILE
C
CEND
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
      DIMENSION QTMP(3)
      LOGICAL TWOCEN, THRCEN, FOUCEN, DERONE, DERTWO
      LOGICAL NODER, DIAGAB, DIAGCD, SHABAB, SHAEQB, SHCEQD,
     *        BIGVEC, DCMPAB, DCMPCD, TCONAB, TCONCD, NOPMAT, NODV
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
#include "mxcent.par"
      DIMENSION ICGTOA(MXAOVC), ICGTOB(MXAOVC),
     *          ICGTOC(MXAOVC), ICGTOD(MXAOVC),IPERT(12)
C
      DIMENSION ICHAIN(NCOORD),NINBCK(NCOORD)
      DIMENSION BUCK(IBKSIZ,NCOORD),IBUCK(IBKSIZ,NCOORD)
      DIMENSION IMAP(NCOORD)
C
      COMMON /INTADR/ IWKAO, IWKSO, IWKHHS, IWK1HH, IWK1HC, IWKLST
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
      DIMENSION SO(1)
C
      COMMON /EXPCOM/  SIGNXYZ(12),
C     *                SIGN1X, SIGN1Y, SIGN1Z, SIGN2X, SIGN2Y, SIGN2Z,
C     *                SIGN3X, SIGN3Y, SIGN3Z, SIGN4X, SIGN4Y, SIGN4Z,
     *                NCENT1, NCENT2, NCENT3, NCENT4,
     *                ISON(4),
c     *                ISO1,   ISO2,   ISO3,   ISO4,
     *                DERONE, DERTWO, TWOCEN, THRCEN, FOUCEN,
     *                NINTYP, NCCINT
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      DIMENSION ITHRPOS(9)
C
      DATA FOURTH,HALF,ONE,TWO,FOUR,ZERADD,AZERO /0.25D0,0.5D0,1.0D0,
     &     2.0D0,4.0D0,1.D-16,0.D0/
      DATA ITHRPOS /1,2,3,1,2,3,1,2,3/
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTXOR(I,J) = XOR(I,J)
C
      ITRI(I,J) = MAX(I,J)*(MAX(I,J) - 1)/2 + MIN(I,J)
      INDEX(I,J,K,L) = NORBA*(NORBB*(NORBC*(L-1)+K-1)+J-1)+I
      IPACK(I,J,K,L)=OR(OR(OR(I,ISHFT(J,IBITWD)),ISHFT(K,2*IBITWD)),
     &                   ISHFT(L,3*IBITWD))

C
C
C  NUMBER OF FIRST DERIVATIVE INTEGRALS
C
      NUMDER=3
      IF(TWOCEN) NUMDER=3
      IF(THRCEN) NUMDER=6
      IF(FOUCEN) NUMDER=9
C
      THRSH=MAX(THRESH,ZERADD)
C
C     ***** Loop over shell components *****
C
      IAOFF=1
      DO 100 ICOMPA=1,KHKTA
       KHKTBB=KHKTB
       IF(DIAGAB) KHKTBB=ICOMPA
       DO 200 ICOMPB = 1,KHKTBB
        DCMPAB=DIAGAB.AND.ICOMPA.EQ.ICOMPB
        DO 300 ICOMPC=1,KHKTC
         KHKTDD=KHKTD
         IF(DIAGCD) KHKTDD=ICOMPC
         DO 400 ICOMPD=1,KHKTDD
          DCMPCD=DIAGCD.AND.ICOMPC.EQ.ICOMPD
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
            DO 700 IORBC=1,NORBC
             INDXC=ICGTOC(IORBC)+IORBC+(ICOMPC-1)*NORBC
             NORBDD=NORBD
             IF(TCONCD) NORBDD=IORBC
             DO 800 IORBD=1,NORBDD
              INDXD=ICGTOD(IORBD)+IORBD+(ICOMPD-1)*NORBD
              INT   = IAOADR
C
C GET INFORMATION ABOUT PERTURBATIONS
C
              IPERT(1)=IMAP(IPTCNT(3*NCENT1-2,0))
              IPERT(2)=IMAP(IPTCNT(3*NCENT1-1,0))
              IPERT(3)=IMAP(IPTCNT(3*NCENT1,0))
              IPERT(4)=IMAP(IPTCNT(3*NCENT2-2,0))
              IPERT(5)=IMAP(IPTCNT(3*NCENT2-1,0))
              IPERT(6)=IMAP(IPTCNT(3*NCENT2,0))
              IF(THRCEN.OR.FOUCEN)THEN
               IPERT(7)=IMAP(IPTCNT(3*NCENT3-2,0))
               IPERT(8)=IMAP(IPTCNT(3*NCENT3-1,0))
               IPERT(9)=IMAP(IPTCNT(3*NCENT3,0))
              ENDIF
              IF(FOUCEN)THEN
               IPERT(10)=IMAP(IPTCNT(3*NCENT4-2,0))
               IPERT(11)=IMAP(IPTCNT(3*NCENT4-1,0))
               IPERT(12)=IMAP(IPTCNT(3*NCENT4,0))
              ENDIF
C
              CALL ZERO(QTMP,3)
             
#ifdef _DEBUG
              if (int+NCCINT*NUMDER.gt.lwork1) then
                 print *, '@DPLUNK: Out-of-bounds read.'
                 call errex
              end if
#endif
              DO 900 I = 1, NUMDER
C
               VALUE=SO(INT)
               QTMP(ITHRPOS(I))=QTMP(ITHRPOS(I))+VALUE
               IF(ABS(VALUE).LT.1.E-09) GO TO 890
C
C PUT THIS INTEGRAL INTO A BUCKET
C
               NBUCK=IPERT(I)
               IF(NBUCK.EQ.0) GO TO 890
               NINBCK(NBUCK)=NINBCK(NBUCK)+1
               BUCK(NINBCK(NBUCK),NBUCK)=VALUE
               IBUCK(NINBCK(NBUCK),NBUCK)=
     &                          IPACK(INDXA,INDXB,INDXC,INDXD)
               IF(NINBCK(NBUCK).EQ.IBKSIZ)THEN
                CALL PLUNK(30,BUCK(1,NBUCK),IBUCK(1,NBUCK),
     &                     ICHAIN(NBUCK),NINBCK(NBUCK),IBKSIZ,
     &                     NREC)
                NINBCK(NBUCK)=0
               ENDIF
C
C NOW DEAL WITH TRANSLATIONAL INVARIANCE
C
    
 
890             CONTINUE
                INT = INT + NCCINT
900            CONTINUE
           DO 901 I=1,3
               IF(ABS(QTMP(I)).LT.1.E-09) GO TO 891
               NBUCK=IPERT(I+NUMDER)
               IF(NBUCK.EQ.0) GO TO 891
               NINBCK(NBUCK)=NINBCK(NBUCK)+1
               BUCK(NINBCK(NBUCK),NBUCK)=-QTMP(I)
               IBUCK(NINBCK(NBUCK),NBUCK)=
     &                          IPACK(INDXA,INDXB,INDXC,INDXD)
               IF(NINBCK(NBUCK).EQ.IBKSIZ)THEN
                CALL PLUNK(30,BUCK(1,NBUCK),IBUCK(1,NBUCK),
     &                     ICHAIN(NBUCK),NINBCK(NBUCK),IBKSIZ,
     &                     NREC)
                NINBCK(NBUCK)=0
               ENDIF
891          CONTINUE
901         CONTINUE
C
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
      RETURN
      END
