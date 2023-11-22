      SUBROUTINE GIINT2(NHKTA,NHKTB,KHKTA,KHKTB,ISTEPA,ISTEPB,
     &                 ISTEPU,ISTEPV,NAHGTF,NATOMC,LDIAG,
     &                 CORPX,CORPY,CORPZ,CORAX,CORAY,CORAZ,
     &                 DIFFABX,DIFFABY,DIFFABZ,
     &                 DIFFPAX,DIFFPAY,DIFFPAZ,TWOG,CHARGE,WORK1)
C
C THIS PROGRAM CALCULATES VARIOUS ONE-ELECTRON INTEGRALS
C REQUIRED FOR THE EVALUATION OF MAGNETIC PROPERTIES.
C THIS ROUTINE DEALS WITH THE TWO-CENTER CASES.
C
CEND
C
C UNIVERSITY OF KARLSRUHE,   AUGUST/91 JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LDIAG
C
#include "baslims.par"
#include "mxcent.par"
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN)
C
      DIMENSION WORK1(1),CHARGE(1)
C
#include "lwork.par"
      PARAMETER (LWKRLM = LWORK3 - 8020)
      COMMON /CWORK3/ WK3LOW, SHGTF, AHGTF((MXCENT+1)*(2*MXQN+1)**3)
CSSS                  RLMCOF(LWKRLM), WK3HGH
#include "odcs.com"
      COMMON /ADER/ ADER0 (MXAQNS)
      COMMON /POINTER/IS0000, IS000X, IS000Y, IS000Z,
     *                IS00XX, IS00XY, IS00XZ, IS00YY,
     *                IS00YZ, IS00ZZ, IT0000, IT000X,
     *                IT000Y, IT000Z, IT00XX, IT00XY,
     *                IT00XZ, IT00YY, IT00YZ, IT00ZZ,
     *                IL0000, IL000X, IL000Y, IL000Z,
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
      COMMON/SCRPOIN/ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,
     *               ISCR9,ISCR10,ISCR11,ISCR12,ISCR13,ISCR14,ISCR15,
     *               ISCR16,ISCR17,ISCR18,ISCR19,ISCR20,
     *               JSCR1,JSCR2,JSCR3,JSCR4,JSCR5,JSCR6,JSCR7,JSCR8,
     *               JSCR9,JSCR10,JSCR11,JSCR12,JSCR13,JSCR14,JSCR15,
     *               JSCR16,JSCR17,JSCR18,JSCR19,JSCR20,JSCR21,JSCR22,
     *               JSCR23,JSCR24,JSCR25,JSCR26,JSCR27,JSCR28 
      COMMON /CENTC/ SIGNCX(MXCENT), SIGNCY(MXCENT), SIGNCZ(MXCENT),
     *               NCENTC(MXCENT), JSYMC(MXCENT),  JCENTC(MXCENT),
     *               ICXVEC(MXCENT), ICYVEC(MXCENT), ICZVEC(MXCENT)
      COMMON /LMNS/ LVALUA(MXAQN), MVALUA(MXAQN), NVALUA(MXAQN),
     *              LVALUB(MXAQN), MVALUB(MXAQN), NVALUB(MXAQN)
      COMMON/GENCON/NRCA,NRCB,CONTA(MXCONT),CONTB(MXCONT)
C
      DATA HALFM/-0.5D0/,TWO /2.D0/,ONEM/-1.D0/,DP5 /0.5D0/ 
      DATA THRSH /1.D-12/
C
C
      DO 100 I = 1,KHKTA*KHKTB
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
C
C CALCULATE DERIVATIVES OF OVERLAP AND ONE-ELECTRON HAMILTONIAN
C INTEGRALS WITH RESPECT TO AN EXTERNAL MAGNETIC FIELD
C
       SX0=SHGTF*ODC00X(ISTRET)
       SY0=SHGTF*ODC00Y(ISTREU)
       SZ0=SHGTF*ODC00Z(ISTREV)
       DX0=SHGTF*ODC00X(ISTRET+1)+CORPX*SX0
       DY0=SHGTF*ODC00Y(ISTREU+1)+CORPY*SY0
       DZ0=SHGTF*ODC00Z(ISTREV+1)+CORPZ*SZ0
       DX00=SHGTF*ODC00X(ISTRET+1)+DIFFPAX*SX0
       DY00=SHGTF*ODC00Y(ISTREU+1)+DIFFPAY*SY0
       DZ00=SHGTF*ODC00Z(ISTREV+1)+DIFFPAZ*SZ0
       SX1=SHGTF*ODC10X(ISTRET)
       SY1=SHGTF*ODC10Y(ISTREU)
       SZ1=SHGTF*ODC10Z(ISTREV)
       SX2=SHGTF*ODC20X(ISTRET)
       SY2=SHGTF*ODC20Y(ISTREU)
       SZ2=SHGTF*ODC20Z(ISTREV)
       DX2=SHGTF*ODC20X(ISTRET+1)+CORPX*SX2
       DY2=SHGTF*ODC20Y(ISTREU+1)+CORPY*SY2
       DZ2=SHGTF*ODC20Z(ISTREV+1)+CORPZ*SZ2
C
       WORK1(ISCR1+I)=SX0*SY0*SZ0
       WORK1(ISCR2+I)=SX0*(-DIFFABY*SY0*DZ0+DIFFABZ*SZ0*DY0)
       WORK1(ISCR3+I)=SY0*(-DIFFABZ*SZ0*DX0+DIFFABX*SX0*DZ0)
       WORK1(ISCR4+I)=SZ0*(-DIFFABX*SX0*DY0+DIFFABY*SY0*DX0)
C
C KE INTEGRALS AND DERIVATIVES OF THE KE INTEGRALS
C
C
C ANGULAR MOMENTUM INTEGRALS
C
C    SIGN: -1 FROM FORMULA, -1 FROM S..1 
C
       WORK1(ISCR5+I)=HALFM*(SX2*SY0*SZ0+SX0*SY2*SZ0+SX0*SY0*SZ2)
C
C  <1|T|0> + <0|T|1> CONTRIBUTION (x component) :
C
       WORK1(ISCR6+I)=HALFM*(-DIFFABY*(DZ0*(SY2*SX0+SY0*SX2)
     &                       +SY0*SX0*DZ2)
     &                     +DIFFABZ*(DY0*(SZ2*SX0+SZ0*SX2) 
     &                       +SZ0*SX0*DY2))
C
C <0|dH/dB|0> CONTRIBUTION (x component) :     
C
     &                 +SX0*(DY00*SZ1-SY1*DZ00)
C
C  <1|T|0> + <0|T|1> CONTRIBUTION (y component) :
C
       WORK1(ISCR7+I)=HALFM*(-DIFFABZ*(DX0*(SY2*SZ0+SY0*SZ2)
     &                       +SY0*SZ0*DX2)
     &                     +DIFFABX*(DZ0*(SX2*SY0+SX0*SY2)
     &                       +SY0*SX0*DZ2))
C
C <0|dH/dB|0> CONTRIBUTION (y component) :    
C
     &                    +SY0*(DZ00*SX1-SZ1*DX00)
C
C  <1|T|0> + <0|T|1> CONTRIBUTION (z component) :
C
       WORK1(ISCR8+I)=HALFM*(-DIFFABX*(DY0*(SZ2*SX0+SZ0*SX2)
     &                       +SZ0*SX0*DY2)
     &                     +DIFFABY*(DX0*(SZ2*SY0+SZ0*SY2)
     &                       +SZ0*SY0*DX2))
C
C <0|dH/dB|0> CONTRIBUTION (z component) :    
C
     &                   +SZ0*(DX00*SY1-SX1*DY00)
C
100   CONTINUE
C
C  CALCULATE INTEGRALS WITH RESPECT TO NUCLEAR MAGNETIC MOMENTS 
C
         MAXADD=2
         INT20 = - NATOMC
         DO 1250 I=1,KHKTA*KHKTB
          ICOMPA=(I-1)/KHKTB+1
          ICOMPB=I-(ICOMPA-1)*KHKTB
          LVALA=LVALUA(ICOMPA)
          MVALA=MVALUA(ICOMPA)
          NVALA=NVALUA(ICOMPA)
          ISTRAT=1+ISTEPA*LVALA
          ISTRAU=1+ISTEPA*MVALA
          ISTRAV=1+ISTEPA*NVALA
          LVALB=LVALUB(ICOMPB)
          MVALB=MVALUB(ICOMPB) 
          NVALB=NVALUB(ICOMPB)
          ISTRET=ISTRAT+ISTEPB*LVALB
          ISTREU=ISTRAU+ISTEPB*MVALB
          ISTREV=ISTRAV+ISTEPB*NVALB
          MAXT = LVALA + LVALB + MAXADD
          MAXU = MVALA + MVALB + MAXADD
          MAXV = NVALA + NVALB + MAXADD
          INT20 = INT20 + NATOMC
          IADRAV = 1
          DO 1200 IV = 0,MAXV
           IADREV = ISTREV + IV
           EV = ODC00Z(IADREV)
           FV = ODC10Z(IADREV)
           IF(IV.EQ.0) THEN
            HV=ODC00Z(IADREV+1)+CORPZ*ODC00Z(IADREV)
            PV=ODC10Z(IADREV+1)+CORPZ*ODC10Z(IADREV)
           ELSE
            HV=(IV+1)*ODC00Z(IADREV+1)+TWOG*ODC00Z(IADREV-1)+
     &         CORPZ*ODC00Z(IADREV)
            PV=(IV+1)*ODC10Z(IADREV+1)+TWOG*ODC10Z(IADREV-1)+
     &         CORPZ*ODC10Z(IADREV)
           ENDIF
           RV=HV-CORAZ*ODC00Z(IADREV)
           IADRAU = IADRAV
           DO 1300 IU = 0,MAXU
            IADREU = ISTREU + IU
            EU = ODC00Y(IADREU)
            FU = ODC10Y(IADREU)
            IF(IU.EQ.0) THEN
             HU=ODC00Y(IADREU+1)+CORPY*ODC00Y(IADREU)
             PU=ODC10Y(IADREU+1)+CORPY*ODC10Y(IADREU)
            ELSE
             HU=(IU+1)*ODC00Y(IADREU+1)+TWOG*ODC00Y(IADREU-1)+
     &          CORPY*ODC00Y(IADREU)
             PU=(IU+1)*ODC10Y(IADREU+1)+TWOG*ODC10Y(IADREU-1)+
     &          CORPY*ODC10Y(IADREU)
            ENDIF
            RU=HU-CORAY*ODC00Y(IADREU)
            EE = EU*EV
            FE = FU*EV
            EF = EU*FV
            FF = FU*FV
            EH = EU*HV
            HE = HU*EV
            FH = FU*HV
            HF = HU*FV
            PE = PU*EV
            EP = EU*PV
            ER = EU*RV
            RE = RU*EV
C
            DO 1400 IT = 0,MAXT
             IADRET = ISTRET + IT
             IADR00 = IADRAU + IT
             IADR0T = IADR00 + 1
             IADR0U = IADR00 + ISTEPU
             IADR0V = IADR00 + ISTEPV
             ET = ODC00X(IADRET)
             FT = ODC10X(IADRET)
             IF(IT.EQ.0) THEN
              HT=ODC00X(IADRET+1)+CORPX*ODC00X(IADRET)
              PT=ODC10X(IADRET+1)+CORPX*ODC10X(IADRET)
             ELSE
              HT=(IT+1)*ODC00X(IADRET+1)+TWOG*ODC00X(IADRET-1)+
     &           CORPX*ODC00X(IADRET)
              PT=(IT+1)*ODC10X(IADRET+1)+TWOG*ODC10X(IADRET-1)+
     &           CORPX*ODC10X(IADRET)
             ENDIF
             RT=HT-CORAX*ODC00X(IADRET)
C
             EEE = ET*EE
             FEE = FT*EE
             EFE = ET*FE
             EEF = ET*EF
             EEH = ET*EH
             EHE = ET*HE
             HEE = HT*EE
             HFE = HT*FE
             HEF = HT*EF
             FHE = FT*HE
             FEH = FT*EH
             EFH = ET*FH
             EHF = ET*HF
             FFE = FT*FE
             FEF = FT*EF
             EFF = ET*FF
             PEE = PT*EE
             EPE = ET*PE
             EEP = ET*EP
             REE = RT*EE
             ERE = ET*RE
             EER = ET*ER
C
             IADRTT = IADR0T + 1
             IADRTU = IADR0T + ISTEPU
             IADRTV = IADR0T + ISTEPV
             IADRUU = IADR0U + ISTEPU
             IADRUV = IADR0U + ISTEPV
             IADRVV = IADR0V + ISTEPV
             IADD = - NAHGTF
C
C                  LOOP OVER NUCLEI 
C
*VOCL LOOP,NOVREC
CDIR$ IVDEP
                  DO 1500 IATOM = 1,NATOMC
C
C                    PICK UP HGTF INTEGRALS
C
                     FACT=ONEM/CHARGE(JCENTC(IATOM))
                     IADD = IADD + NAHGTF
                     AH00 = AHGTF(IADR00 + IADD)
                     AH0T = AHGTF(IADR0T + IADD)
                     AH0U = AHGTF(IADR0U + IADD)
                     AH0V = AHGTF(IADR0V + IADD)
C
C    DERIVATIVES OF NUCLEAR ATTRACTION INTEGRAL
C    WITH RESPECT TO AN EXTERNAL MAGNETIC FIELD :
C
C     <1|V|0> + <0|V|1> CONTRIBUTIONS (x,y,z components) :
C
C  BX
                     WORK1(ISCR6+I)=WORK1(ISCR6+I)+
     &                  AH00*(DIFFABZ*EHE-DIFFABY*EEH)
C  BY
                     WORK1(ISCR7+I)=WORK1(ISCR7+I)+
     &                  AH00*(DIFFABX*EEH-DIFFABZ*HEE)
C  BZ
                     WORK1(ISCR8+I)=WORK1(ISCR8+I)+
     &                  AH00*(DIFFABY*HEE-DIFFABX*EHE)
C
                     INT2 = INT20 + IATOM
C
C DERIVATIVES WITH RESPECT TO THE NUCLEAR MAGNETIC MOMENTS
C
C  INX
                     WORK1(JSCR1+INT2)=WORK1(JSCR1+INT2)+FACT*
     &                                 (AH0U*EEF-AH0V*EFE)
C  INY
                     WORK1(JSCR2+INT2)=WORK1(JSCR2+INT2)+FACT*
     &                                 (AH0V*FEE-AH0T*EEF)
C  INZ
                     WORK1(JSCR3+INT2)=WORK1(JSCR3+INT2)+FACT*
     &                                 (AH0T*EFE-AH0U*FEE)
C
C SECOND DERIVATIVES WITH RESPECT TO NUCLEAR MAGNETIC MOMENTS
C AND EXTERNAL MAGNETIC FIELD
C FIRST INDEX REFERS TO EXTERNAL FIELD
C
C   BX AND INX
                     WORK1(JSCR4+INT2)=WORK1(JSCR4+INT2)+FACT*
     &                                 (-DIFFABY*AH0U*EEP+
     &                                  DIFFABY*AH0V*EFH+
     &                                  DIFFABZ*AH0U*EHF-
     &                                  DIFFABZ*AH0V*EPE-
     &                                  AH0U*ERE-AH0V*EER)
C   BX AND INY
                     WORK1(JSCR5+INT2)=WORK1(JSCR5+INT2)+FACT*
     &                                 (-DIFFABY*AH0V*FEH+
     &                                  DIFFABY*AH0T*EEP+
     &                                  DIFFABZ*AH0V*FHE-
     &                                  DIFFABZ*AH0T*EHF+
     &                                  AH0T*ERE)
C   BX AND INZ
                     WORK1(JSCR6+INT2)=WORK1(JSCR6+INT2)+FACT*
     &                                 (-DIFFABY*AH0T*EFH+
     &                                  DIFFABY*AH0U*FEH+
     &                                  DIFFABZ*AH0T*EPE-
     &                                  DIFFABZ*AH0U*FHE+
     &                                  AH0T*EER)
C   BY AND INX
                     WORK1(JSCR7+INT2)=WORK1(JSCR7+INT2)+FACT*
     &                                 (-DIFFABZ*AH0U*HEF+
     &                                  DIFFABZ*AH0V*HFE+
     &                                  DIFFABX*AH0U*EEP-
     &                                  DIFFABX*AH0V*EFH+
     &                                  AH0U*REE)
C   BY AND INY
                     WORK1(JSCR8+INT2)=WORK1(JSCR8+INT2)+FACT*
     &                                 (-DIFFABZ*AH0V*PEE+
     &                                  DIFFABZ*AH0T*HEF+
     &                                  DIFFABX*AH0V*FEH-
     &                                  DIFFABX*AH0T*EEP-
     &                                  AH0T*REE-AH0V*EER)
C   BY AND INZ
                     WORK1(JSCR9+INT2)=WORK1(JSCR9+INT2)+FACT*
     &                                 (-DIFFABZ*AH0T*HFE+
     &                                  DIFFABZ*AH0U*PEE+
     &                                  DIFFABX*AH0T*EFH-
     &                                  DIFFABX*AH0U*FEH+
     &                                  AH0U*EER)
C   BZ AND INX
                     WORK1(JSCR10+INT2)=WORK1(JSCR10+INT2)+FACT*
     &                                 (-DIFFABX*AH0U*EHF+
     &                                  DIFFABX*AH0V*EPE+
     &                                  DIFFABY*AH0U*HEF-
     &                                  DIFFABY*AH0V*HFE+
     &                                  AH0V*REE)
C   BZ AND INY
                     WORK1(JSCR11+INT2)=WORK1(JSCR11+INT2)+FACT*
     &                                 (-DIFFABX*AH0V*FHE+
     &                                  DIFFABX*AH0T*EHF+
     &                                  DIFFABY*AH0V*PEE-
     &                                  DIFFABY*AH0T*HEF+
     &                                  AH0V*ERE)
C   BZ AND INZ
                     WORK1(JSCR12+INT2)=WORK1(JSCR12+INT2)+FACT*
     &                                 (-DIFFABX*AH0T*EPE+
     &                                  DIFFABX*AH0U*FHE+
     &                                  DIFFABY*AH0T*HFE-
     &                                  DIFFABY*AH0U*PEE-
     &                                  AH0T*REE-AH0U*ERE)
C
 1500             CONTINUE
 1400          CONTINUE
               IADRAU = IADRAU + ISTEPU
 1300       CONTINUE
            IADRAV = IADRAV + ISTEPV
 1200    CONTINUE
 1250 CONTINUE
C
C
C  NOW FORM THE CONTRACTIONS
C
      ISKIP=NRCA*NRCB
      IF(LDIAG) ISKIP=NRCA*(NRCA+1)/2
      IOFF=0
      DO 150 IRCA=1,NRCA
       MAXB=NRCB
       IF(LDIAG) MAXB=IRCA
       DO 150 IRCB=1,MAXB
        IOFF=IOFF+1
        CONT=CONTA(IRCA)*CONTB(IRCB)
        IF(ABS(CONT).GT.THRSH) THEN
         CALL SAXPY(KHKTA*KHKTB*8,CONT,WORK1(ISCR1+1),1,
     &              WORK1(IS0000+IOFF),ISKIP)
        ENDIF
150   CONTINUE
      INUMB=3
      ISKIP=NRCA*NRCB*NATOMC
      IF(LDIAG) ISKIP=NRCA*(NRCA+1)/2*NATOMC
      IOFF=0
      DO 1700 IRCA=1,NRCA
       MAXB=NRCB
       IF(LDIAG) MAXB=IRCA
       DO 1700 IRCB=1,MAXB
        CONT=CONTA(IRCA)*CONTB(IRCB)
        IF(ABS(CONT).GT.THRSH) THEN
         DO 1710 IATOM=1,NATOMC
          CALL SAXPY(KHKTA*KHKTB*INUMB,CONT,WORK1(JSCR1+IATOM),NATOMC,
     &               WORK1(IA000X+IOFF+IATOM),ISKIP)
1710     CONTINUE
        ENDIF
        IOFF=IOFF+NATOMC
1700  CONTINUE
C
      INUMB=9
      IOFF=0
      DO 1701 IRCA=1,NRCA
       MAXB=NRCB
       IF(LDIAG) MAXB=IRCA
       DO 1701 IRCB=1,MAXB
        CONT=CONTA(IRCA)*CONTB(IRCB)
        IF(ABS(CONT).GT.THRSH) THEN
         DO 1711 IATOM=1,NATOMC
          CALL SAXPY(KHKTA*KHKTB*INUMB,CONT,WORK1(JSCR4+IATOM),NATOMC,
     &               WORK1(IA0X0X+IOFF+IATOM),ISKIP)
1711     CONTINUE
        ENDIF
        IOFF=IOFF+NATOMC
1701  CONTINUE
C
      RETURN
      END
