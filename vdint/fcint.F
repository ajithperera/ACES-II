      SUBROUTINE FCINT(NHKTA,NHKTB,KHKTA,KHKTB,LDIAG,
     &                 XA,YA,ZA,EXPA,XB,YB,ZB,EXPB,
     &                 CORCX,CORCY,CORCZ,FACINT,
     &                 NATOMC,WORK1)
C
C THIS ROUTINE CALCULATED THE FERMI-CONTACT INTEGRALS
C 
CEND
C
C 4/93 JG KARLSRUHE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LDIAG
C
#include "baslims.par"
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN
     &           *MXCONT*MXCONT)
C
      DIMENSION WORK1(1)
      DIMENSION CORCX(NATOMC),CORCY(NATOMC),CORCZ(NATOMC)
      DIMENSION FACINT(NATOMC)
C
      COMMON /POINTER/IS0000, IS000X, IS000Y, IS000Z,
     &                IS00XX, IS00XY, IS00XZ, IS00YY,
     &                IS00YZ, IS00ZZ, IT0000, IT000X,
     &                IT000Y, IT000Z, IT00XX, IT00XY,
     &                IT00XZ, IT00YY, IT00YZ, IT00ZZ,
     &                ID0000, ID000X, ID000Y, ID000Z,
     &                ID00XX, ID00XY, ID00XZ, ID00YX,
     &                ID00YY, ID00YZ, ID00ZX, ID00ZY,
     &                ID00ZZ,
     &                IA0000, IA0X00, IA0Y00, IA0Z00,
     &                IAXX00, IAXY00, IAXZ00, IAYY00,
     &                IAYZ00, IAZZ00, IA000X, IA000Y,
     &                IA000Z, IA00XX, IA00XY, IA00XZ,
     &                IA00YY, IA00YZ, IA00ZZ, IA0X0X,
     &                IA0X0Y, IA0X0Z, IA0Y0X, IA0Y0Y,
     &                IA0Y0Z, IA0Z0X, IA0Z0Y, IA0Z0Z
      COMMON/SCRPOIN/ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,
     &               ISCR9,ISCR10,ISCR11,ISCR12,ISCR13,ISCR14,ISCR15,
     &               ISCR16,ISCR17,ISCR18,ISCR19,ISCR20,
     &               JSCR1,JSCR2,JSCR3,JSCR4,JSCR5,JSCR6,JSCR7,JSCR8,
     &               JSCR9,JSCR10,JSCR11,JSCR12,JSCR13,JSCR14,JSCR15,
     &               JSCR16,JSCR17,JSCR18,JSCR19,JSCR20,JSCR21,JSCR22,
     &               JSCR23,JSCR24,JSCR25,JSCR26,JSCR27,JSCR28 
#include "mxcent.par"
      COMMON /CENTC/ SIGNCX(MXCENT), SIGNCY(MXCENT), SIGNCZ(MXCENT),
     &               NCENTC(MXCENT), JSYMC(MXCENT),  JCENTC(MXCENT),
     &               ICXVEC(MXCENT), ICYVEC(MXCENT), ICZVEC(MXCENT)
      COMMON /LMNS/ LVALUA(MXAQN), MVALUA(MXAQN), NVALUA(MXAQN),
     &              LVALUB(MXAQN), MVALUB(MXAQN), NVALUB(MXAQN)
      COMMON/GENCON/NRCA,NRCB,CONTA(MXCONT),CONTB(MXCONT)
C
      DATA THRSH /1.D-20/
C
C CALCULATE NUCLEAR ATTRACTION INTEGRALS *****
C
      INT20 = - NATOMC
      DO 250 I=1,KHKTA*KHKTB
       ICOMPA=(I-1)/KHKTB+1
       ICOMPB=I-(ICOMPA-1)*KHKTB
       LVALA=LVALUA(ICOMPA)
       MVALA=MVALUA(ICOMPA)
       NVALA=NVALUA(ICOMPA)
       LVALB=LVALUB(ICOMPB)
       MVALB=MVALUB(ICOMPB) 
       NVALB=NVALUB(ICOMPB)
       INT20 = INT20 + NATOMC
C
C LOOP OVER NUCLEI 
C
*VOCL LOOP,NOVREC
CDIR$ IVDEP
       DO 500 IATOM = 1,NATOMC
C
        FACT=FACINT(IATOM)
        XC=CORCX(IATOM)
        YC=CORCY(IATOM)
        ZC=CORCZ(IATOM)
        RAC=(XA-XC)**2+(YA-YC)**2+(ZA-ZC)**2
        RBC=(XB-XC)**2+(YB-YC)**2+(ZB-ZC)**2
C
        INT2 = INT20 + IATOM
C
C FC INTEGRALS
C
        if(lvala.eq.0)then
         fact1=1.0d0
        else
         fact1=(xc-xa)**lvala
        endif 

        if(lvalb.eq.0)then
         fact2=1.0d0
        else
         fact2=(xc-xb)**lvalb
        endif 

        if(mvala.eq.0)then
         fact3=1.0d0
        else
         fact3=(yc-ya)**mvala
        endif 

        if(mvalb.eq.0)then
         fact4=1.0d0
        else
         fact4=(yc-yb)**mvalb
        endif 

        if(nvala.eq.0)then
         fact5=1.0d0
        else
         fact5=(zc-za)**nvala
        endif 

        if(nvalb.eq.0)then
         fact6=1.0d0
        else
         fact6=(zc-zb)**nvalb
        endif 

        WORK1(JSCR1+INT2)= fact1*fact2*fact3*fact4*fact5*fact6
     &                    *EXP(-EXPA*RAC-EXPB*RBC)*FACT
C 
500    CONTINUE
250   CONTINUE
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
     &              WORK1(IS000X+IOFF+IATOM),ISKIP)
710     CONTINUE
        ENDIF
       IOFF=IOFF+NATOMC
700   CONTINUE
C
C ALL DONE, RETURN
      RETURN
      END
