      SUBROUTINE CINT1(NHKTA,NHKTB,KHKTA,KHKTB,ICENTA,LDIAG,
     *                 ISTEPA,ISTEPB,ISTEPU,ISTEPV,NAHGTF,
     *                 NATOMC,SECOND,WORK1)
C
C     TUH
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LDIAG, SECOND
      DIMENSION WORK1(1)
#include "lwork.par"
      PARAMETER (LWKRLM = LWORK3 - 8020)
#include "baslims.par"
#include "mxcent.par"
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN
     *           *MXCONT*MXCONT)
#include "odcs.com"
      COMMON /CWORK3/ WK3LOW, SHGTF, AHGTF((MXCENT+1)*(2*MXQN+1)**3)
CSSS                  RLMCOF(LWKRLM),WK3HGH
#include "odcs.com"
      COMMON /ADER/ ADER0 (MXAQNS)
      COMMON /POINTER/IS0000, IS000X, IS000Y, IS000Z,
     *                IS00XX, IS00XY, IS00XZ, IS00YY,
     *                IS00YZ, IS00ZZ, IT0000, IT000X,
     *                IT000Y, IT000Z, IT00XX, IT00XY,
     *                IT00XZ, IT00YY, IT00YZ, IT00ZZ,
     *                ID0000, ID000X, ID000Y, ID000Z,
     *                ID00XX, ID00XY, ID00XZ, ID00YX,
     *                ID00YY, ID00YZ, ID00ZX, ID00ZY,
     *                ID00ZZ,
     *                IA0000, IA0X00, IA0Y00, IA0Z00,
     *                IAXX00, IAXY00, IAXZ00, IAYY00,
     *                IAYZ00, IAZZ00, IA000X, IA000Y,
     *                IA000Z, IA00XX, IA00XY, IA00XZ,
     *                IA00YY, IA00YZ, IA00ZZ, IA0X0X,
     *                IA0X0Y, IA0X0Z, IA0Y0X, IA0Y0Y,
     *                IA0Y0Z, IA0Z0X, IA0Z0Y, IA0Z0Z
      COMMON/SCRPOIN/ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,
     *               ISCR8,ISCR9,ISCR10,ISCR11,ISCR12,ISCR13,ISCR14,
     *               ISCR15,ISCR16,ISCR17,ISCR18,ISCR19,ISCR20,
     *               JSCR1,JSCR2,JSCR3,JSCR4,JSCR5,JSCR6,
     *               JSCR7,JSCR8,JSCR9,JSCR10,JSCR11,JSCR12,JSCR13,
     *               JSCR14,JSCR15,JSCR16,JSCR17,JSCR18,JSCR19,
     *               JSCR20,JSCR21,JSCR22,JSCR23,JSCR24,JSCR25,JSCR26,
     *               JSCR27,JSCR28 
      COMMON /CENTC/ SIGNCX(MXCENT), SIGNCY(MXCENT), SIGNCZ(MXCENT),
     *               NCENTC(MXCENT), JSYMC(MXCENT),  JCENTC(MXCENT),
     *               ICXVEC(MXCENT), ICYVEC(MXCENT), ICZVEC(MXCENT)
      COMMON /LMNS/ LVALUA(MXAQN), MVALUA(MXAQN), NVALUA(MXAQN),
     *              LVALUB(MXAQN), MVALUB(MXAQN), NVALUB(MXAQN)
C
      COMMON /GENCON/ NRCA,NRCB,CONTA(MXCONT),CONTB(MXCONT)
C
      DATA THRSH /1.D-20/
C
      INT20=-NATOMC
      DO 100 I=1,KHKTA*KHKTB
       ICOMPA=(I-1)/KHKTB+1
       ICOMPB=I-(ICOMPA-1)*KHKTB
       LVALA = LVALUA(ICOMPA)
       MVALA = MVALUA(ICOMPA)
       NVALA = NVALUA(ICOMPA)
       ISTRAT = 1 + ISTEPA*LVALA
       ISTRAU = 1 + ISTEPA*MVALA
       ISTRAV = 1 + ISTEPA*NVALA
       LVALB = LVALUB(ICOMPB)
       MVALB = MVALUB(ICOMPB)
       NVALB = NVALUB(ICOMPB)
       ISTRET = ISTRAT + ISTEPB*LVALB
       ISTREU = ISTRAU + ISTEPB*MVALB
       ISTREV = ISTRAV + ISTEPB*NVALB
       MAXT = LVALA + LVALB
       MAXU = MVALA + MVALB
       MAXV = NVALA + NVALB
       INT20 = INT20 + NATOMC
       IADRAV = 1
       DO 200 IV = 0,MAXV
        EV = ODC00Z(ISTREV + IV)
        IADRAU = IADRAV
        DO 300 IU = 0,MAXU
         EE = ODC00Y(ISTREU + IU)*EV
         DO 400 IT = 0,MAXT
          EEE = ODC00X(ISTRET + IT)*EE
          IADR00 = IADRAU + IT
          IADR0T = IADR00 + 1
          IADR0U = IADR00 + ISTEPU
          IADR0V = IADR00 + ISTEPV
          IF (SECOND) THEN
           IADRTT = IADR0T + 1
           IADRTU = IADR0T + ISTEPU
           IADRTV = IADR0T + ISTEPV
           IADRUU = IADR0U + ISTEPU
           IADRUV = IADR0U + ISTEPV
           IADRVV = IADR0V + ISTEPV
          END IF
          IADD = - NAHGTF
C
C   LOOP OVER NUCLEI 
C
*VOCL LOOP,NOVREC
CDIR$ IVDEP
          DO 500 IATOM = 1,NATOMC
           IADD = IADD + NAHGTF
           ICENTC = NCENTC(IATOM)
           IF (ICENTC .NE. ICENTA) THEN
            INT2 = INT20 + IATOM
            AH00 = AHGTF(IADR00 + IADD)
            AH0T = AHGTF(IADR0T + IADD)
            AH0U = AHGTF(IADR0U + IADD)
            AH0V = AHGTF(IADR0V + IADD)
C
C  UNDIFFERENTIATED INTEGRAL:
C
            WORK1(JSCR1+INT2)=WORK1(JSCR1+INT2)+EEE*AH00
c                        WORK1(IA0000+INT2) = WORK1(IA0000+INT2)+EEE*AH00
C
C  C DIFFERENTIATED INTEGRALS:
C
            WORK1(JSCR2+INT2)=WORK1(JSCR2+INT2)-EEE*AH0T
            WORK1(JSCR3+INT2)=WORK1(JSCR3+INT2)-EEE*AH0U
            WORK1(JSCR4+INT2)=WORK1(JSCR4+INT2)-EEE*AH0V

C            Write(6,"(3(1x,F20.10))") WORK1(JSCR2+INT2),
C     &                                WORK1(JSCR3+INT2),
C     &                                WORK1(JSCR4+INT2)

c                        WORK1(IA000X+INT2) = WORK1(IA000X+INT2)-EEE*AH0T
c                        WORK1(IA000Y+INT2) = WORK1(IA000Y+INT2)-EEE*AH0U
c                        WORK1(IA000Z+INT2) = WORK1(IA000Z+INT2)-EEE*AH0V
C
C                       SECOND DERIVATIVES:
C
            IF (SECOND) THEN
                        AHTT = AHGTF(IADRTT + IADD)
                           AHTU = AHGTF(IADRTU + IADD)
                           AHTV = AHGTF(IADRTV + IADD)
                           AHUU = AHGTF(IADRUU + IADD)
                           AHUV = AHGTF(IADRUV + IADD)
                           AHVV = AHGTF(IADRVV + IADD)
C
C                          C-C DIFFERENTIATED INTEGRALS:
C
c                           WORK1(IA00XX + INT2) = WORK1(IA00XX + INT2)
c     *                                          + EEE*AHTT
c                           WORK1(IA00XY + INT2) = WORK1(IA00XY + INT2)
c     *                                          + EEE*AHTU
c                           WORK1(IA00XZ + INT2) = WORK1(IA00XZ + INT2)
c     *                                          + EEE*AHTV
c                           WORK1(IA00YY + INT2) = WORK1(IA00YY + INT2)
c     *                                          + EEE*AHUU
c                           WORK1(IA00YZ + INT2) = WORK1(IA00YZ + INT2)
c     *                                          + EEE*AHUV
c                           WORK1(IA00ZZ + INT2) = WORK1(IA00ZZ + INT2)
c     *                                          + EEE*AHVV
                           WORK1(JSCR5 + INT2) = WORK1(JSCR5 + INT2)
     *                                          + EEE*AHTT
                           WORK1(JSCR6 + INT2) = WORK1(JSCR6 + INT2)
     *                                          + EEE*AHTU
                           WORK1(JSCR7 + INT2) = WORK1(JSCR7 + INT2)
     *                                          + EEE*AHTV
                           WORK1(JSCR8+ INT2) = WORK1(JSCR8 + INT2)
     *                                          + EEE*AHUU
                           WORK1(JSCR9 + INT2) = WORK1(JSCR9 + INT2)
     *                                          + EEE*AHUV
                           WORK1(JSCR10 + INT2) = WORK1(JSCR10 + INT2)
     *                                          + EEE*AHVV
                        END IF
                     END IF
  500             CONTINUE
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
  100 CONTINUE
C
C  NOW FORM THE CONTRACTIONS
C
      ISKIP=NRCA*NRCB*NATOMC
      IF(LDIAG) ISKIP=NRCA*(NRCA+1)/2*NATOMC
      IOFF=0
      DO 700 IRCA=1,NRCA
       MAXB=NRCB
       IF(LDIAG) MAXB=IRCA
       DO 700 IRCB=1,MAXB
        CONT=CONTA(IRCA)*CONTB(IRCB)
        IF(ABS(CONT).GT.THRSH) THEN
         DO 710 IATOM=1,NATOMC
         CALL SAXPY(KHKTA*KHKTB,CONT,WORK1(JSCR1+IATOM),NATOMC,
     *              WORK1(IA0000+IOFF+IATOM),ISKIP)
         CALL SAXPY(KHKTA*KHKTB,CONT,WORK1(JSCR2+IATOM),NATOMC,
     *             WORK1(IA000X+IOFF+IATOM),ISKIP)
         CALL SAXPY(KHKTA*KHKTB,CONT,WORK1(JSCR3+IATOM),NATOMC,
     *             WORK1(IA000Y+IOFF+IATOM),ISKIP)
         CALL SAXPY(KHKTA*KHKTB,CONT,WORK1(JSCR4+IATOM),NATOMC,
     *             WORK1(IA000Z+IOFF+IATOM),ISKIP)
         IF(SECOND) THEN
         CALL SAXPY(KHKTA*KHKTB,CONT,WORK1(JSCR5+IATOM),NATOMC,
     *             WORK1(IA00XX+IOFF+IATOM),ISKIP)
         CALL SAXPY(KHKTA*KHKTB,CONT,WORK1(JSCR6+IATOM),NATOMC,
     *             WORK1(IA00XY+IOFF+IATOM),ISKIP)
         CALL SAXPY(KHKTA*KHKTB,CONT,WORK1(JSCR7+IATOM),NATOMC,
     *             WORK1(IA00XZ+IOFF+IATOM),ISKIP)
         CALL SAXPY(KHKTA*KHKTB,CONT,WORK1(JSCR8+IATOM),NATOMC,
     *             WORK1(IA00YY+IOFF+IATOM),ISKIP)
         CALL SAXPY(KHKTA*KHKTB,CONT,WORK1(JSCR9+IATOM),NATOMC,
     *             WORK1(IA00YZ+IOFF+IATOM),ISKIP)
         CALL SAXPY(KHKTA*KHKTB,CONT,WORK1(JSCR10+IATOM),NATOMC,
     *             WORK1(IA00ZZ+IOFF+IATOM),ISKIP)
         ENDIF
710     CONTINUE
        ENDIF
       IOFF=IOFF+NATOMC
700   CONTINUE

      RETURN
      END
