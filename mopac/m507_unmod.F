      FUNCTION AABABC(IOCCA1, IOCCB1, IOCCA2, IOCCB2, NMOS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION IOCCA1(NMOS), IOCCB1(NMOS), IOCCA2(NMOS), IOCCB2(NMOS)
***********************************************************************
*
* AABABC EVALUATES THE C.I. MATRIX ELEMENT FOR TWO MICROSTATES DIFFERING
*       BY BETA ELECTRON. THAT IS, ONE MICROSTATE HAS A BETA ELECTRON
*       IN PSI(I) WHICH, IN THE OTHER MICROSTATE IS IN PSI(J)
*
***********************************************************************
      COMMON /XYIJKL/ XY(NMECI,NMECI,NMECI,NMECI)
      COMMON /BASEOC/ OCCA(NMECI)
         SAVE                                                           GL0892
      DO 10 I=1,NMOS
   10 IF(IOCCA1(I).NE.IOCCA2(I)) GOTO 20
   20 IJ=IOCCB1(I)
      DO 30 J=I+1,NMOS
         IF(IOCCA1(J).NE.IOCCA2(J)) GOTO 40
   30 IJ=IJ+IOCCA1(J)+IOCCB1(J)
   40 SUM=0.D0
      DO 50 K=1,NMOS
   50 SUM=SUM+ (XY(I,J,K,K)-XY(I,K,J,K))*(IOCCA1(K)-OCCA(K)) +
     1          XY(I,J,K,K)             *(IOCCB1(K)-OCCA(K))
      AABABC=SUM*((-1)**(IJ-(IJ/2)*2))
      RETURN
      END
      FUNCTION AABACD(IOCCA1, IOCCB1, IOCCA2, IOCCB2, NMOS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION IOCCA1(NMOS), IOCCB1(NMOS), IOCCA2(NMOS), IOCCB2(NMOS)
***********************************************************************
*
* AABACD EVALUATES THE C.I. MATRIX ELEMENT FOR TWO MICROSTATES DIFFERING
*       BY TWO ALPHA MOS. ONE MICROSTATE HAS ALPHA ELECTRONS IN
*       M.O.S PSI(I) AND PSI(J) FOR WHICH THE OTHER MICROSTATE HAS
*       ELECTRONS IN PSI(K) AND PSI(L)
*
***********************************************************************
      COMMON /XYIJKL/ XY(NMECI,NMECI,NMECI,NMECI)
         SAVE                                                           GL0892
      IJ=0
      DO 10 I=1,NMOS
   10 IF(IOCCA1(I) .LT. IOCCA2(I)) GOTO 20
   20 DO 30 J=I+1,NMOS
         IF(IOCCA1(J) .LT. IOCCA2(J)) GOTO 40
   30 IJ=IJ+IOCCA2(J)+IOCCB2(J)
   40 DO 50 K=1,NMOS
   50 IF(IOCCA1(K) .GT. IOCCA2(K)) GOTO 60
   60 DO 70 L=K+1,NMOS
         IF(IOCCA1(L) .GT. IOCCA2(L)) GOTO 80
   70 IJ=IJ+IOCCA1(L)+IOCCB1(L)
   80 IJ=IJ+IOCCB2(I)+IOCCB1(K)
      AABACD=(XY(I,K,J,L)-XY(I,L,K,J))*((-1)**(IJ-(IJ/2)*2))
      RETURN
      END
      FUNCTION AABBCD(IOCCA1, IOCCB1, IOCCA2, IOCCB2, NMOS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION IOCCA1(NMOS), IOCCB1(NMOS), IOCCA2(NMOS), IOCCB2(NMOS)
***********************************************************************
*
* AABBCD EVALUATES THE C.I. MATRIX ELEMENT FOR TWO MICROSTATES DIFFERING
*       BY TWO SETS OF M.O.S. ONE MICROSTATE HAS AN ALPHA ELECTRON
*       IN PSI(I) AND A BETA ELECTRON IN PSI(K) FOR WHICH THE OTHER
*       MICROSTATE HAS AN ALPHA ELECTRON IN PSI(J) AND A BETA ELECTRON
*       IN PSI(L)
*
***********************************************************************
      COMMON /XYIJKL/ XY(NMECI,NMECI,NMECI,NMECI)
      COMMON /SPQR/ ISPQR(NMECI*NMECI,NMECI),IS,ILOOP, JLOOP
         SAVE                                                           GL0892
      DO 10 I=1,NMOS
   10 IF(IOCCA1(I) .NE. IOCCA2(I)) GOTO 20
   20 DO 30 J=I+1,NMOS
   30 IF(IOCCA1(J) .NE. IOCCA2(J)) GOTO 40
   40 DO 50 K=1,NMOS
   50 IF(IOCCB1(K) .NE. IOCCB2(K)) GOTO 60
   60 DO 70 L=K+1,NMOS
   70 IF(IOCCB1(L) .NE. IOCCB2(L)) GOTO 80
   80 IF( I.EQ.K .AND. J.EQ.L .AND. IOCCA1(I).NE.IOCCB1(I)) THEN
         ISPQR(ILOOP,IS)=JLOOP
         IS=IS+1
      ENDIF
      IF(IOCCA1(I) .LT. IOCCA2(I)) THEN
         M=I
         I=J
         J=M
      ENDIF
      IF(IOCCB1(K) .LT. IOCCB2(K)) THEN
         M=K
         K=L
         L=M
      ENDIF
      XR=XY(I,J,K,L)
C#      WRITE(6,'(4I5,F12.6)')I,J,K,L,XR
C
C   NOW UNTANGLE THE MICROSTATES
C
      IJ=1
      IF( I.GT.K .AND. J.GT.L .OR. I.LE.K .AND. J.LE.L)IJ=0
      IF( I.GT.K ) IJ=IJ+IOCCA1(K)+IOCCB1(I)
      IF( J.GT.L ) IJ=IJ+IOCCA2(L)+IOCCB2(J)
      IF(I.GT.K)THEN
         M=I
         I=K
         K=M
      ENDIF
      DO 90 M=I,K
   90 IJ=IJ+IOCCB1(M)+IOCCA1(M)
      IF(J.GT.L)THEN
         M=J
         J=L
         L=M
      ENDIF
      DO 100 M=J,L
  100 IJ=IJ+IOCCB2(M)+IOCCA2(M)
C
C   IJ IN THE PERMUTATION NUMBER, .EQUIV. -1 IF IJ IS ODD.
C
      AABBCD=XR*((-1)**(IJ-(IJ/2)*2))
      RETURN
      END
C
      SUBROUTINE ANAVIB(COORD,EIGS,DIPT,N3,VIBS,RIJ,HESS,TRAVEL,REDMAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION COORD(3,NUMATM),EIGS(N3),VIBS(N3,N3), DIPT(N3),
     1RIJ(MAXHES), TRAVEL(N3), HESS(*), REDMAS(*)
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /ELEMTS/ ELEMNT(107)
      LOGICAL VIB1, VIB2, VIB3, VIB4, VIB5
      COMMON /FOKMAT/ F(MPACK*2)
      CHARACTER*2 ELEMNT
      DIMENSION VANRAD(107), IJF(10), FIJ(10)
         SAVE                                                           GL0892
      DATA VANRAD/
     1   0.32,0.93,
     2   1.23, 0.90, 0.82, 0.77, 0.75, 0.73, 0.72, 0.71,
     3   1.54, 1.36, 1.18, 1.11, 1.06, 1.02, 0.99, 0.98,
     4   2.03, 1.74, 1.44, 1.32, 1.22, 1.18, 1.17, 1.17, 1.16,
     5   1.15, 1.17, 1.25, 1.26, 1.22, 1.20, 1.16, 1.14, 1.12,
     6   2.16, 1.91, 1.62, 1.45, 1.34, 1.30, 1.27, 1.25, 1.25,
     7   1.28, 1.34, 1.48, 1.44, 1.41, 1.40, 1.36, 1.33, 1.31,
     8   2.35, 1.98, 1.69,
     9   1.65, 1.65, 1.64, 1.63, 1.62, 1.85, 1.61, 1.59, 1.59, 1.58,
     1   1.57, 1.56, 1.56, 1.56,
     2   1.44, 1.34, 1.30, 1.28, 1.26, 1.27, 1.30, 1.34,
     3   1.49, 1.48, 1.47, 1.46, 1.46, 1.45,1.45,21*1.45/
      DATA FACT/6.023D23/
      C2PI=1.D0/(2.998D10*3.141592653598D0*2.D0)
      N3=NUMAT*3
C
C    COMPUTE INTERATOMIC DISTANCES.
C
      L=0
      DO 10 I=1,NUMAT
         DO 10 J=1,I-1
            L=L+1
   10 RIJ(L)=SQRT((COORD(1,J)-COORD(1,I))**2+
     1            (COORD(2,J)-COORD(2,I))**2+
     2            (COORD(3,J)-COORD(3,I))**2)+1.D-10
C
C     ANALYSE VIBRATIONS
C
      WRITE(6,'(''1'',//10X,''DESCRIPTION OF VIBRATIONS'',/)')
      ILINE=0
      DO 120 K=1,N3
         IF(ABS(EIGS(K)).LT.50) GOTO 120
         VIB1=.TRUE.
         VIB2=.TRUE.
         VIB3=.TRUE.
         VIB4=.TRUE.
         VIB5=.TRUE.
         J3=0
         L=0
         TOT=0.D0
         LINEAR=0
         J1=-2
         DO 70 J=1,NUMAT
            J1=J1+3
            I1=-2
            DO 70 I=1,J-1
               I1=I1+3
               VDW=(VANRAD(NAT(I))+VANRAD(NAT(J)))*1.5
               L=L+1
               F(L)=0.D0
               IF(   RIJ(L)  .LT.  VDW) THEN
C
C CALCULATE ENERGY TERM BETWEEN THE TWO ATOMS
C
                  EAB=0.D0
                  DO 20 JJ=J1,J1+2
                     DO 20 II=I1,I1+2
   20             EAB=EAB+VIBS(JJ,K)*HESS((JJ*(JJ-1))/2+II)*VIBS(II,K)
                  EB=0.D0
                  DO 40 JJ=J1,J1+2
                     DO 30 II=J1,JJ
   30                EB=EB+VIBS(JJ,K)*HESS((JJ*(JJ-1))/2+II)*VIBS(II,K)*
     12.D0
   40             EB=EB-VIBS(JJ,K)*HESS((JJ*(JJ+1))/2)*VIBS(JJ,K)
                  EA=0.D0
                  DO 60 JJ=I1,I1+2
                     DO 50 II=I1,JJ
   50                EA=EA+VIBS(JJ,K)*HESS((JJ*(JJ-1))/2+II)*VIBS(II,K)*
     12.D0
   60             EA=EA-VIBS(JJ,K)*HESS((JJ*(JJ+1))/2)*VIBS(JJ,K)
                  LINEAR=LINEAR+1
                  F(L)=EA+EAB*2.D0+EB
                  TOT=TOT+F(L)
               ENDIF
   70    CONTINUE
C
C  NOW TO SORT F INTO DECENDING ORDER
C
         DO 90 I=1,10
            SUM=-100.D0
            DO 80 J=1,L
               IF(ABS(F(J)).GT.SUM)THEN
                  JJ=J
                  SUM=ABS(F(J))
               ENDIF
   80       CONTINUE
            IF(SUM.LT.0.D0)GOTO 100
            FIJ(I)=SUM
            F(JJ)=-1.D-5
            IJF(I)=JJ
C#      WRITE(6,*)FIJ(I),IJF(I)
   90    CONTINUE
         I=10
  100    LINEAR=I
         SUM=1.D0/(TOT+1.D-8)
         DO 110 IJ=1,LINEAR
            J=0.5D0*(0.99D0+SQRT(1.D0+8.D0*IJF(IJ)))
            I=IJF(IJ)-(J*(J-1))/2
            J=J+1
            XJ=COORD(1,J)
            YJ=COORD(2,J)
            ZJ=COORD(3,J)
            J1=3*J-2
            J2=J1+1
            J3=J2+1
            I3=0
            XI=COORD(1,I)
            YI=COORD(2,I)
            ZI=COORD(3,I)
            I1=3*I-2
            I2=I1+1
            I3=I2+1
            X= VIBS(J1,K)-VIBS(I1,K)
            Y= VIBS(J2,K)-VIBS(I2,K)
            Z= VIBS(J3,K)-VIBS(I3,K)
            E=FIJ(IJ)*SUM*100.D0
            IF(ABS(E).GT.110)GOTO 120
            SHIFT=X*X+Y*Y+Z*Z+1.D-30
            IF(ABS(E) .GT. 10.D0.OR.IJ.LT.5.AND.ABS(E).GT.0.1) THEN
               SHIFT=SQRT(SHIFT)
               RADIAL=((X*(XI-XJ)+Y*(YI-YJ)+Z*(ZI-ZJ))
     1                  /(SHIFT*RIJ(IJF(IJ))))**2*100.D0
               IF (VIB1) THEN
                  WRITE(6,'(/,'' VIBRATION'',I4,''            ATOM PAIR
     1 '',''    ENERGY CONTRIBUTION              RADIAL'')')K
                  ANS=100.D0*SQRT(FIJ(IJ)*1.D5*6.023D23)/(2.998D10*3.141
     159D0*2.D0)/EIGS(K)
                  ANS=MIN(999.9D0,MAX(-99.9D0,ANS))
                  WRITE(6,'('' FREQ.   '',F9.2,6X,1A2,I2,
     1'' -- '',A2,I2, ''         '',F6.1,''% ('',F5.1,''%)'',F18.1
     2,''%'')')
     3EIGS(K),ELEMNT(NAT(I)),I,ELEMNT(NAT(J)),J,E,ANS,RADIAL
C#      WRITE(6,*)ANS
                  VIB1=.FALSE.
               ELSEIF (VIB2) THEN
                  VIB2=.FALSE.
                  WRITE(6,'('' T-DIPOLE'',F9.4,6X,1A2,I2,
     1'' -- '',A2,I2, ''         '',F6.1,''%'',F27.1,''%'')')
     2DIPT(K),ELEMNT(NAT(I)),I,ELEMNT(NAT(J)),J,E,RADIAL
               ELSEIF (VIB3) THEN
                  VIB3=.FALSE.
                  WRITE(6,'('' TRAVEL  '',F9.4,6X,1A2,I2,
     1'' -- '',A2,I2, ''         '',F6.1,''%'',F27.1,''%'')')
     2TRAVEL(K),ELEMNT(NAT(I)),I,ELEMNT(NAT(J)),J,E,RADIAL
               ELSEIF (VIB4) THEN
                  VIB4=.FALSE.
                  WRITE(6,'('' RED. MASS'',F8.4,6X,1A2,I2,
     1'' -- '',A2,I2, ''         '',F6.1,''%'',F27.1,''%'')')
     2REDMAS(K),ELEMNT(NAT(I)),I,ELEMNT(NAT(J)),J,E,RADIAL
               ELSE
                  ILINE=ILINE+1
                  WRITE(6,'(''                        '',1A2,I2,
     1'' -- '',A2,I2, ''         '',F6.1,''%'',F27.1,''%'')')
     2ELEMNT(NAT(I)),I,ELEMNT(NAT(J)),J,E,RADIAL
               ENDIF
            ENDIF
  110    CONTINUE
         ILINE=ILINE+6
         IF(VIB1)WRITE(6,'(/,'' VIBRATION'',I4)')K
         IF(VIB1)WRITE(6,'(  '' FREQ.    '',F8.2)')EIGS(K)
         IF(VIB2)WRITE(6,'(  '' T-DIPOLE '',F8.4)')DIPT(K)
         IF(VIB3)WRITE(6,'(  '' TRAVEL   '',F8.4)')TRAVEL(K)
         IF(VIB4)WRITE(6,'(  '' RED. MASS'',F8.4)')REDMAS(K)
         IF(ILINE.GT.52)THEN
            ILINE=0
            WRITE(6,'(''1'')')
         ENDIF
  120 CONTINUE
      RETURN
      END
      SUBROUTINE AXIS(COORD,NUMAT,A,B,C,SUMW, MASS,EVEC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION COORD(3,NUMAT)
      COMMON /KEYWRD/ KEYWRD
************************************************************************
*
*  AXIS CALCULATES THE THREE MOMENTS OF INERTIA AND THE MOLECULAR
*       WEIGHT.  THE MOMENTS OF INERTIA ARE RETURNED IN A, B, AND C.
*       THE MOLECULAR WEIGHT IN SUMW.
*       THE UNITS OF INERTIA ARE 10**(-40)GRAM-CM**2,
*       AND MOL.WEIGHT IN ATOMIC-MASS-UNITS. (AMU'S)
************************************************************************
      COMMON /ATMASS/ ATMASS(NUMATM)
      DIMENSION T(6), X(NUMATM), Y(NUMATM),
     1          Z(NUMATM), ROT(3), XYZMOM(3), EIG(3), EVEC(3,3)
      COMMON /NUMCAL/ NUMCAL                                            IR0295
      CHARACTER*80 KEYWRD
         SAVE                                                           GL0892
      DATA T /6*0.D0/, ICALCN /0/                                       IR0295
************************************************************************
*     CONST1 =  10**40/(N*A*A)
*               N = AVERGADRO'S NUMBER
*               A = CM IN AN ANGSTROM
*               10**40 IS TO ALLOW UNITS TO BE 10**(-40)GRAM-CM**2
*
************************************************************************
      CONST1 = 1.66053D0
************************************************************************
*
*     CONST2 = CONVERSION FACTOR FROM ANGSTROM-AMU TO CM**(-1)
*
*            = (PLANCK'S CONSTANT*N*10**16)/(8*PI*PI*C)
*            = 6.62618*10**(-27)[ERG-SEC]*6.02205*10**23*10**16/
*              (8*(3.1415926535)**2*2.997925*10**10[CM/SEC])
*
************************************************************************
      CONST2=16.8576522D0
C    FIRST WE CENTRE THE MOLECULE ABOUT THE CENTRE OF GRAVITY,
C    THIS DEPENDS ON THE ISOTOPIC MASSES, AND THE CARTESIAN GEOMETRY.
C
      SUMW=1.D-20
      SUMWX=0.D0
      SUMWY=0.D0
      SUMWZ=0.D0
      WEIGHT=1.D0
      DO 10 I=1,NUMAT
         IF(MASS.GT.0)WEIGHT=ATMASS(I)
         SUMW=SUMW+WEIGHT
         SUMWX=SUMWX+WEIGHT*COORD(1,I)
         SUMWY=SUMWY+WEIGHT*COORD(2,I)
   10 SUMWZ=SUMWZ+WEIGHT*COORD(3,I)
      SUMWX=SUMWX/SUMW
      SUMWY=SUMWY/SUMW
      SUMWZ=SUMWZ/SUMW
      DO 20 I=1,NUMAT
         X(I)=COORD(1,I)-SUMWX
         Y(I)=COORD(2,I)-SUMWY
   20 Z(I)=COORD(3,I)-SUMWZ
************************************************************************
*
*    MATRIX FOR MOMENTS OF INERTIA IS OF FORM
*
*           |   Y**2+Z**2                         |
*           |    -Y*X       Z**2+X**2             | -I =0
*           |    -Z*X        -Z*Y       X**2+Y**2 |
*
************************************************************************
      DO 30 I=1,6
   30 T(I)=I*1.D-10
      DO 40 I=1,NUMAT
         IF(MASS.GT.0)WEIGHT=ATMASS(I)
         T(1)=T(1)+WEIGHT*(Y(I)**2+Z(I)**2)
         T(2)=T(2)-WEIGHT*X(I)*Y(I)
         T(3)=T(3)+WEIGHT*(Z(I)**2+X(I)**2)
         T(4)=T(4)-WEIGHT*Z(I)*X(I)
         T(5)=T(5)-WEIGHT*Y(I)*Z(I)
   40 T(6)=T(6)+WEIGHT*(X(I)**2+Y(I)**2)
      CALL RSP(T,3,3,EIG,EVEC)
      IF(MASS.GT.0.AND. (ICALCN .NE. NUMCAL).AND.
     &   INDEX(KEYWRD,'RC=').EQ.0) THEN
         WRITE(6,'(/10X,''MOLECULAR WEIGHT ='',F8.2,/)')
     &             MIN(99999.99D0,SUMW)
         WRITE(6,'(//10X,'' PRINCIPAL MOMENTS OF INERTIA IN CM(-1)'',/)'
     1)
         DO 50 I=1,3
            IF(EIG(I).LT.3.D-4) THEN
               EIG(I)=0.D0
               ROT(I)=0.D0
            ELSE
               ROT(I)=CONST2/EIG(I)
            ENDIF
   50    XYZMOM(I)=EIG(I)*CONST1
         WRITE(6,'(10X,''A ='',F12.6,''   B ='',F12.6,
     1''   C ='',F12.6,/)')(ROT(I),I=1,3)
         IF(INDEX(KEYWRD,'RC=').EQ.0)
     1WRITE(6,'(//10X,'' PRINCIPAL MOMENTS OF INERTIA IN '',
     2''UNITS OF 10**(-40)*GRAM-CM**2'',/)')
         WRITE(6,'(10X,''A ='',F12.6,''   B ='',F12.6,
     1''   C ='',F12.6,/)')(XYZMOM(I),I=1,3)
         C=ROT(1)
         B=ROT(2)
         A=ROT(3)
      ENDIF
C
C   NOW TO ORIENT THE MOLECULE SO THE CHIRALITY IS PRESERVED
C
      SUM=EVEC(1,1)*(EVEC(2,2)*EVEC(3,3)-EVEC(3,2)*EVEC(2,3)) +
     1    EVEC(1,2)*(EVEC(2,3)*EVEC(3,1)-EVEC(2,1)*EVEC(3,3)) +
     2    EVEC(1,3)*(EVEC(2,1)*EVEC(3,2)-EVEC(2,2)*EVEC(3,1))
      IF( SUM .LT. 0) THEN
         DO 60 J=1,3
   60    EVEC(J,1)=-EVEC(J,1)
      ENDIF
      DO 70 I=1,NUMAT
         COORD(1,I)=X(I)
         COORD(2,I)=Y(I)
         COORD(3,I)=Z(I)
   70 CONTINUE
      IF(MASS.GT.0) ICALCN=NUMCAL                                        IR0295
      END
      FUNCTION BABBBC(IOCCA1, IOCCB1, IOCCA2, IOCCB2, NMOS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
      DIMENSION IOCCA1(NMOS), IOCCB1(NMOS), IOCCA2(NMOS), IOCCB2(NMOS)
***********************************************************************
*
* BABBBC EVALUATES THE C.I. MATRIX ELEMENT FOR TWO MICROSTATES DIFFERING
*       BY ONE BETA ELECTRON. THAT IS, ONE MICROSTATE HAS A BETA
*       ELECTRON IN PSI(I) AND THE OTHER MICROSTATE HAS AN ELECTRON IN
*       PSI(J).
***********************************************************************
      COMMON /XYIJKL/ XY(NMECI,NMECI,NMECI,NMECI)
      COMMON /BASEOC/ OCCA(NMECI)
         SAVE                                                           GL0892
      DO 10 I=1,NMOS
   10 IF(IOCCB1(I).NE.IOCCB2(I)) GOTO 20
   20 IJ=0
      DO 30 J=I+1,NMOS
         IF(IOCCB1(J).NE.IOCCB2(J)) GOTO 40
   30 IJ=IJ+IOCCA1(J)+IOCCB1(J)
   40 IJ=IJ+IOCCA1(J)
C
C   THE UNPAIRED M.O.S ARE I AND J
      SUM=0.D0
      DO 50 K=1,NMOS
   50 SUM=SUM+ (XY(I,J,K,K)-XY(I,K,J,K))*(IOCCB1(K)-OCCA(K)) +
     1          XY(I,J,K,K)             *(IOCCA1(K)-OCCA(K))
      BABBBC=SUM*((-1)**(IJ-(IJ/2)*2))
      RETURN
      END
      FUNCTION BABBCD(IOCCA1, IOCCB1, IOCCA2, IOCCB2, NMOS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION IOCCA1(NMOS), IOCCB1(NMOS), IOCCA2(NMOS), IOCCB2(NMOS)
***********************************************************************
*
* BABBCD EVALUATES THE C.I. MATRIX ELEMENT FOR TWO MICROSTATES DIFFERING
*       BY TWO BETA MOS. ONE MICROSTATE HAS BETA ELECTRONS IN
*       M.O.S PSI(I) AND PSI(J) FOR WHICH THE OTHER MICROSTATE HAS
*       ELECTRONS IN PSI(K) AND PSI(L)
*
***********************************************************************
      COMMON /XYIJKL/ XY(NMECI,NMECI,NMECI,NMECI)
         SAVE                                                           GL0892
      IJ=0
      DO 10 I=1,NMOS
   10 IF(IOCCB1(I) .LT. IOCCB2(I)) GOTO 20
   20 DO 30 J=I+1,NMOS
         IF(IOCCB1(J) .LT. IOCCB2(J)) GOTO 40
   30 IJ=IJ+IOCCA2(J)+IOCCB2(J)
   40 IJ=IJ+IOCCA2(J)
      DO 50 K=1,NMOS
   50 IF(IOCCB1(K) .GT. IOCCB2(K)) GOTO 60
   60 DO 70 L=K+1,NMOS
         IF(IOCCB1(L) .GT. IOCCB2(L)) GOTO 80
   70 IJ=IJ+IOCCA1(L)+IOCCB1(L)
   80 IJ=IJ+IOCCA1(L)
      IF((IJ/2)*2.EQ.IJ) THEN
         ONE=1.D0
      ELSE
         ONE=-1.D0
      ENDIF
      BABBCD=(XY(I,K,J,L)-XY(I,L,J,K))*ONE
      RETURN
      END
      SUBROUTINE BANGLE(XYZ,I,J,K,ANGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XYZ(3,*)
         SAVE                                                           GL0892
*********************************************************************
*
* BANGLE CALCULATES THE ANGLE BETWEEN ATOMS I,J, AND K. THE
*        CARTESIAN COORDINATES ARE IN XYZ.
*
*********************************************************************
      D2IJ = (XYZ(1,I)-XYZ(1,J))**2+
     1       (XYZ(2,I)-XYZ(2,J))**2+
     2       (XYZ(3,I)-XYZ(3,J))**2
      D2JK = (XYZ(1,J)-XYZ(1,K))**2+
     1       (XYZ(2,J)-XYZ(2,K))**2+
     2       (XYZ(3,J)-XYZ(3,K))**2
      D2IK = (XYZ(1,I)-XYZ(1,K))**2+
     1       (XYZ(2,I)-XYZ(2,K))**2+
     2       (XYZ(3,I)-XYZ(3,K))**2
      XY = SQRT(D2IJ*D2JK)
      TEMP = 0.5D0 * (D2IJ+D2JK-D2IK) / XY
      IF (TEMP .GT. 1.0D0) TEMP=1.0D0
      IF (TEMP .LT. -1.0D0) TEMP=-1.0D0
      ANGLE = ACOS( TEMP )
      RETURN
      END
      SUBROUTINE BFN(X,BF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BF(13)
C**********************************************************************
C
C     BINTGS FORMS THE "B" INTEGRALS FOR THE OVERLAP CALCULATION.
C
C**********************************************************************
      DIMENSION FACT(17)
         SAVE                                                           GL0892
      DATA FACT/1.D0,2.D0,6.D0,24.D0,120.D0,720.D0,5040.D0,40320.D0,
     1362880.D0,3628800.D0,39916800.D0,479001600.D0,6227020800.D0,
     28.71782912D10,1.307674368D12,2.092278989D13,3.556874281D14/
      K=12
      IO=0
      ABSX = ABS(X)
      IF (ABSX.GT.3.D00) GO TO 40
      IF (ABSX.LE.2.D00) GO TO 10
      LAST=15
      GO TO 60
   10 IF (ABSX.LE.1.D00) GO TO 20
      LAST=12
      GO TO 60
   20 IF (ABSX.LE.0.5D00) GO TO 30
      LAST=7
      GO TO 60
   30 IF (ABSX.LE.1.D-6) GOTO 90
      LAST=6
      GO TO 60
   40 EXPX=EXP(X)
      EXPMX=1.D00/EXPX
      xinv = 1.0/x
      bf(1) = (expx - expmx)*xinv
      zk = 1.0
      DO 50 I=1,K
      zk = -1.0*zk
      bf(i+1) = (i*bf(i) + zk*expx - expmx)*xinv
50    continue
      return
   60 DO 80 I=IO,K
          zk = -x
          y = 2.0*mod(i+1,2)/(i+1)
          do 70 m=io+1,last
          y = y + 2.0*zk*mod(m+i+1,2)/(fact(m)*(m+i+1))
          zk = -x*zk
70        continue
   80 BF(I+1)=Y
      return
   90 DO 100 I=IO,K
  100 BF(I+1)=(2*MOD(I+1,2))/(I+1.D0)
      RETURN
C
      END
      SUBROUTINE CALPAR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /ONELEC/ USS(107),UPP(107),UDD(107)
     1       /ATOMIC/ EISOL(107),EHEAT(107)
     2       /ALPHA / ALP(107)
     3       /EXPONT/ ZS(107),ZP(107),ZD(107)
     4       /GAUSS / FN1(107),FN2(107)
     5       /BETAS / BETAS(107),BETAP(107),BETAD(107)
     6       /TWOELE/ GSS(107),GSP(107),GPP(107),GP2(107),HSP(107),
     7                GSD(107),GPD(107),GDD(107)
     8       /IDEAS / GUESS1(107,10), GUESS2(107,10), GUESS3(107,10)
      COMMON /MNDO/  USSM(107), UPPM(107), UDDM(107), ZSM(107),ZPM(107),
     1ZDM(107), BETASM(107), BETAPM(107), BETADM(107), ALPM(107),
     2EISOLM(107), DDM(107), QQM(107), AMM(107), ADM(107), AQM(107)
     3,GSSM(107),GSPM(107),GPPM(107),GP2M(107),HSPM(107), POLVOM(107)
      COMMON /KEYWRD/ KEYWRD
      COMMON/MULTIP/ DD(107),QQ(107),AM(107),AD(107),AQ(107)
      DIMENSION NSPQN(107)
      CHARACTER KEYWRD*80
      DIMENSION USSC(107), UPPC(107), GSSC(107), GSPC(107), HSPC(107),
     1GP2C(107), GPPC(107), UDDC(107), GSDC(107), GDDC(107)
         SAVE                                                           GL0892
      DATA NSPQN/2*1,8*2,8*3,18*4,18*5,32*6,21*0/
C
C THE CONTINUATION LINES INDICATE THE PRINCIPAL QUANTUM NUMBER.
C
      DATA USSC/
     11.D0,                                                      0.D0,
     21.D0,                                               6*2.D0,0.D0,
     31.D0,                                               6*2.D0,0.D0,
     41.D0,4*2.D0,1.D0,4*2.D0,1.D0,                       6*2.D0,0.D0,
     51.D0,3*2.D0,2*1.D0,2.D0,2*1.D0,0.D0,1.D0,           6*2.D0,0.D0,
     61.D0,22*2.D0,1.D0,1.D0,                             6*2.D0,0.D0,
     721*0.D0/
      DATA  UPPC/
     1 2*0.D0,
     2 2*0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,
     3 2*0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,
     412*0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,
     512*0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,
     626*0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,
     721*0.D0/
      DATA UDDC/18*0.D0,
     1 2*0.D0,1.D0,2.D0,3.D0,5.D0,5.D0,6.D0,7.D0,8.D0,1.D1,1.D1,6*0.D0,
     2 2*0.D0,1.D0,2.D0,4.D0,5.D0,5.D0,7.D0,8.D0,1.D1,1.D1,1.D1,6*0.D0,
     3 2*0.D0,1.D0,6*0.D0,1.D0,6*0.D0,1.D0,2.D0,3.D0,4.D0,
     4                           5.D0,6.D0,7.D0,9.D0,1.D1,1.D1,6*0.D0,
     5 21*0.D0/
      DATA GSSC/2*0.D0,
     1 0.D0,6*1.D0,0.D0,
     2 0.D0,6*1.D0,0.D0,
     3 0.D0,4*1.D0,0.D0,4*1.D0,0.D0,6*1.D0,0.D0,
     4 0.D0,3*1.D0,7*0.D0,6*1.D0,0.D0,
     5 0.D0,22*1.D0,2*1.D0,6*1.D0,0.D0,
     6 21*0.D0/
      DATA GSPC/2*0.D0,
     1 2*0.D0,2.D0,4.D0,6.D0,8.D0,10.D0,0.D0,
     2 2*0.D0,2.D0,4.D0,6.D0,8.D0,10.D0,0.D0,
     312*0.D0,2.D0,4.D0,6.D0,8.D0,10.D0,0.D0,
     412*0.D0,2.D0,4.D0,6.D0,8.D0,10.D0,0.D0,
     526*0.D0,2.D0,4.D0,6.D0,8.D0,10.D0,0.D0,
     621*0.D0/
      DATA HSPC/2*0.D0,
     1 2*0.D0,-1.D0,-2.D0,-3.D0,-4.D0,-5.D0,0.D0,
     2 2*0.D0,-1.D0,-2.D0,-3.D0,-4.D0,-5.D0,0.D0,
     312*0.D0,-1.D0,-2.D0,-3.D0,-4.D0,-5.D0,0.D0,
     412*0.D0,-1.D0,-2.D0,-3.D0,-4.D0,-5.D0,0.D0,
     526*0.D0,-1.D0,-2.D0,-3.D0,-4.D0,-5.D0,0.D0,
     621*0.D0/
      DATA GP2C/2*0.D0,
     1 3*0.D0,1.5D0,4.5D0,6.5D0,10.D0,0.D0,
     2 3*0.D0,1.5D0,4.5D0,6.5D0,10.D0,0.D0,
     313*0.D0,1.5D0,4.5D0,6.5D0,10.D0,0.D0,
     413*0.D0,1.5D0,4.5D0,6.5D0,10.D0,0.D0,
     527*0.D0,1.5D0,4.5D0,6.5D0,10.D0,0.D0,
     621*0.D0/
      DATA GPPC/2*0.D0,
     1 3*0.D0,-0.5D0,-1.5D0,-0.5D0,2*0.D0,
     2 3*0.D0,-0.5D0,-1.5D0,-0.5D0,2*0.D0,
     313*0.D0,-0.5D0,-1.5D0,-0.5D0,2*0.D0,
     413*0.D0,-0.5D0,-1.5D0,-0.5D0,2*0.D0,
     527*0.D0,-0.5D0,-1.5D0,-0.5D0,2*0.D0,
     621*0.D0/
     7GSDC/18*0.D0,
     8 2*0.D0,2.D0,4.D0,6.D0,5.D0,10.D0,12.D0,14.D0,16.D0,10.D0,7*0.D0,
     9 2*0.D0,2.D0,4.D0,4.D0,5.D0,6.D0,7.D0,8.D0,0.D0,1.D1,7*0.D0,
     1 2*0.D0,2.D0,6*0.D0,2.D0,6*0.D0,2.D0,4.D0,6.D0,8.D0,10.D0,12.D0,
     2                                         14.D0,9.D0,10.D0,7*0.D0,
     321*0.D0/
      DATA GDDC/18*0.D0,
     1 3*0.D0,1.D0,3.D0,10.D0,10.D0,15.D0,21.D0,28.D0,8*0.D0,
     2 3*0.D0,1.D0,6.D0,10.D0,15.D0,21.D0,28.D0,45.D0,8*0.D0,
     317*0.D0,1.D0,3.D0, 6.D0,10.D0,15.D0,21.D0,36.D0,8*0.D0,
     421*0.D0/
C
C  The DATA block shown above is derived from the ground-state atomic
C  configuration of the elements.  In checking it, pay careful attention
C  to the actual ground-state configuration. Note also that there are no
C  configurations which have both p and d electrons in the valence shell
C
C     SET SCALING PARAMETER.
      P=2.D0
      P2=P*P
      P4=P**4
      DO 30 I=2,98
         IF(ZP(I).LT.1.D-4.AND.ZS(I).LT.1.D-4)GOTO 30
**********************************************************************
*
*   CONSTRAINTS ON THE POSSIBLE VALUES OF PARAMETERS
*
**********************************************************************
         IF(ZP(I).LT.0.3D0) ZP(I)=0.3D0
C  PUT IN ANY CONSTRAINTS AT THIS POINT
**********************************************************************
         HPP=0.5D0*(GPP(I)-GP2(I))
         HPP=MAX(0.1D0,HPP)
         HSP(I)=MAX(1.D-7,HSP(I))
         EISOL(I)=USS(I)*USSC(I)+UPP(I)*UPPC(I)+UDD(I)*UDDC(I)+
     1         GSS(I)*GSSC(I)+GPP(I)*GPPC(I)+GSP(I)*GSPC(I)+
     2         GP2(I)*GP2C(I)+HSP(I)*HSPC(I)+GSD(I)*GSDC(I)+
     3         GDD(I)*GDDC(I)
         QN=NSPQN(I)
         DD(I)=(2.D0*QN+1)*(4.D0*ZS(I)*ZP(I))**(QN+0.5D0)/(ZS(I)+ZP(I))
     1**(2.D0*QN+2)/SQRT(3.D0)
         DDM(I)=DD(I)
         QQ(I)=SQRT((4.D0*QN*QN+6.D0*QN+2.D0)/20.D0)/ZP(I)
         QQM(I)=QQ(I)
C     CALCULATE ADDITIVE TERMS, IN ATOMIC UNITS.
         JMAX=5
         GDD1= (P2*HSP(I)/(27.21* 4.*DD(I)**2))**(1./3.)
         GQQ= (P4*HPP/(27.21*48.*QQ(I)**4))**0.2
         D1=GDD1
         D2=GDD1+0.04
         Q1=GQQ
         Q2=GQQ+0.04
         DO 10 J=1,JMAX
            DF=D2-D1
            HSP1= 2.*D1 - 2./SQRT(4.*DD(I)**2+1./D1**2)
            HSP2= 2.*D2 - 2./SQRT(4.*DD(I)**2+1./D2**2)
            HSP1= HSP1/P2
            HSP2= HSP2/P2
            D3= D1 + DF*(HSP(I)/27.21-HSP1)/(HSP2-HSP1)
            D1= D2
            D2= D3
   10    CONTINUE
         DO 20 J=1,JMAX
            QF=Q2-Q1
            HPP1= 4.*Q1 - 8./SQRT(4.*QQ(I)**2+1./Q1**2)
     1            + 4./SQRT(8.*QQ(I)**2+1./Q1**2)
            HPP2= 4.*Q2 - 8./SQRT(4.*QQ(I)**2+1./Q2**2)
     1            + 4./SQRT(8.*QQ(I)**2+1./Q2**2)
            HPP1= HPP1/P4
            HPP2= HPP2/P4
            Q3= Q1 + QF*(HPP/27.21-HPP1)/(HPP2-HPP1)
            Q1= Q2
            Q2= Q3
   20    CONTINUE
         AM(I)= GSS(I)/27.21
         AD(I)= D2
         AQ(I)= Q2
         AMM(I)=AM(I)
         ADM(I)=AD(I)
         AQM(I)=AQ(I)
   30 CONTINUE
      EISOL(1)=USS(1)
      AM(1)=GSS(1)/27.21D0
      AD(1)=AM(1)
      AQ(1)=AM(1)
      AMM(1)=AM(1)
      ADM(1)=AD(1)
      AQM(1)=AQ(1)
C
C     DEBUG PRINTING.
C     THIS IS FORMATTED FOR DIRECT INSERTION INTO 'PARAM'
C
      IF(INDEX(KEYWRD,'DEP').EQ.0) RETURN
      WRITE(6,50)
      DO 60 I=1,107
         IF(ZS(I).EQ.0) GOTO 60
         WRITE(6,'(''C'',20X,''DATA FOR ELEMENT'',I3)')I
         WRITE(6,'(6X,''DATA USS   ('',I3,'')/'',F16.7,''D0/'')')
     1                    I,USS(I)
         IF(UPP(I) .NE. 0.D0)
     1WRITE(6,'(6X,''DATA UPP   ('',I3,'')/'',F16.7,''D0/'')')I,UPP(I)
         IF(UDD(I) .NE. 0.D0)
     1WRITE(6,'(6X,''DATA UDD   ('',I3,'')/'',F16.7,''D0/'')')I,UDD(I)
         IF(BETAS(I) .NE. 0.D0)
     1WRITE(6,'(6X,''DATA BETAS ('',I3,'')/'',F16.7,''D0/'')')
     2I,BETAS(I)
         IF(BETAP(I) .NE. 0.D0)
     1WRITE(6,'(6X,''DATA BETAP ('',I3,'')/'',F16.7,''D0/'')')
     2I,BETAP(I)
         IF(BETAD(I) .NE. 0.D0)
     1WRITE(6,'(6X,''DATA BETAD ('',I3,'')/'',F16.7,''D0/'')')
     2I,BETAD(I)
         WRITE(6,'(6X,''DATA ZS    ('',I3,'')/'',F16.7,''D0/'')')
     1I,ZS(I)
         IF(ZP(I) .NE. 0.D0)
     1WRITE(6,'(6X,''DATA ZP    ('',I3,'')/'',F16.7,''D0/'')')I,ZP(I)
         IF(ZD(I) .NE. 0.D0)
     1WRITE(6,'(6X,''DATA ZD    ('',I3,'')/'',F16.7,''D0/'')')I,ZD(I)
         WRITE(6,'(6X,''DATA ALP   ('',I3,'')/'',F16.7,''D0/'')')
     1I,ALP(I)
         WRITE(6,'(6X,''DATA EISOL ('',I3,'')/'',F16.7,''D0/'')')
     1I,EISOL(I)
         IF(GSS(I) .NE. 0.D0)
     1WRITE(6,'(6X,''DATA GSS   ('',I3,'')/'',F16.7,''D0/'')')
     2I,GSS(I)
         IF(GSP(I) .NE. 0.D0)
     1WRITE(6,'(6X,''DATA GSP   ('',I3,'')/'',F16.7,''D0/'')')
     2I,GSP(I)
         IF(GPP(I) .NE. 0.D0)
     1WRITE(6,'(6X,''DATA GPP   ('',I3,'')/'',F16.7,''D0/'')')
     2I,GPP(I)
         IF(GP2(I) .NE. 0.D0)
     1WRITE(6,'(6X,''DATA GP2   ('',I3,'')/'',F16.7,''D0/'')')
     2I,GP2(I)
         IF(HSP(I) .NE. 0.D0)
     1WRITE(6,'(6X,''DATA HSP   ('',I3,'')/'',F16.7,''D0/'')')
     2I,HSP(I)
         IF(DD(I) .NE. 0.D0)
     1WRITE(6,'(6X,''DATA DD    ('',I3,'')/'',F16.7,''D0/'')')I,DD(I)
         IF(QQ(I) .NE. 0.D0)
     1WRITE(6,'(6X,''DATA QQ    ('',I3,'')/'',F16.7,''D0/'')')I,QQ(I)
         WRITE(6,'(6X,''DATA AM    ('',I3,'')/'',F16.7,''D0/'')')
     1I,AM(I)
         IF(AD(I) .NE. 0.D0)
     1WRITE(6,'(6X,''DATA AD    ('',I3,'')/'',F16.7,''D0/'')')I,AD(I)
         IF(AQ(I) .NE. 0.D0)
     1WRITE(6,'(6X,''DATA AQ    ('',I3,'')/'',F16.7,''D0/'')')I,AQ(I)
         IF(FN1(I) .NE. 0.D0)
     1WRITE(6,'(6X,''DATA FN1   ('',I3,'')/'',F16.7,''D0/'')')I,FN1(I)
         IF(FN2(I) .NE. 0.D0)
     1WRITE(6,'(6X,''DATA FN2   ('',I3,'')/'',F16.7,''D0/'')')I,FN2(I)
         DO 40 J=1,10
            IF(GUESS1(I,J) .EQ.0.D0) GOTO 40
            WRITE(6,'(6X,''DATA GUESS1('',I3,'','',I1,'')/'',
     1F16.7,''D0/'')')I,J,GUESS1(I,J)
            WRITE(6,'(6X,''DATA GUESS2('',I3,'','',I1,'')/'',
     1F16.7,''D0/'')')I,J,GUESS2(I,J)
            WRITE(6,'(6X,''DATA GUESS3('',I3,'','',I1,'')/'',
     1F16.7,''D0/'')')I,J,GUESS3(I,J)
   40    CONTINUE
   50    FORMAT(1H ,1X,'OUTPUT INCLUDES DEBUG INFORMATION',//)
   60 CONTINUE
      RETURN
      END
      FUNCTION CAPCOR(NAT,NFIRST,NLAST,NUMAT,P,H)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION P(*), H(*), NFIRST(NUMAT), NLAST(NUMAT), NAT(NUMAT)
         SAVE                                                           GL0892
******************************************************************
*
*    CORRECTION TO ELECTRONIC ENERGY DUE TO CAPPED BONDS
*
******************************************************************
      SUM=0.D0
      DO 40 I=1,NUMAT
         NI=NAT(I)
         IL=NFIRST(I)
         IU=NLAST(I)
         IF(NI.EQ.102) THEN
C
C   DO ENTIRE ROW - NO NEED TO CHECK FURTHER.
C
            J=(NLAST(I)*(NLAST(I)+1))/2
            II=IU-1
            DO 10 K=1,II
               J=J-1
   10       SUM=SUM+P(J)*H(J)
         ELSE
            DO 30 J=1,I
               JL=NFIRST(J)
               IF(NAT(J).EQ.102)THEN
                  DO  20 K=IL,IU
                     KK=(K*(K-1))/2+JL
   20             SUM=SUM+P(KK)*H(KK)
               ENDIF
   30       CONTINUE
         ENDIF
   40 CONTINUE
C
C   DOUBLE SUM SINCE WE ONLY CALCULATED LOWER HALF, AND CAPCOR
C   WILL APPEAR IN 1/2*P(H+F).  ONLY H PART OF F WILL BE USED.
      CAPCOR=-SUM*2.D0
      RETURN
      END
      SUBROUTINE CHRGE(P,Q)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION P(*),Q(*)
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
         SAVE                                                           GL0892
C***********************************************************************
C
C      CHRGE STORES IN Q THE TOTAL ELECTRON DENSITIES ON THE ATOMS
C
C      ON INPUT P      = DENSITY MATRIX
C
C      ON OUTPUT Q     = ATOM ELECTRON DENSITIES
C
C***********************************************************************
      K=0
      DO 10 I=1,NUMAT
         IA=NFIRST(I)
         IB=NLAST(I)
         Q(I)=0.D0
         DO 10 J=IA,IB
            K=K+J
   10 Q(I)=Q(I)+P(K)
      RETURN
      END
      SUBROUTINE COE(X1,Y1,Z1,X2,Y2,Z2,PQ1,PQ2,C,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER PQ1,PQ2,PQ,CO
      DIMENSION C(75)
         SAVE                                                           GL0892
      XY=(X2-X1)**2+(Y2-Y1)**2
      R=SQRT(XY+(Z2-Z1)**2)
      XY=SQRT(XY)
      IF (XY.LT.1.D-10) GO TO 10
      CA=(X2-X1)/XY
      CB=(Z2-Z1)/R
      SA=(Y2-Y1)/XY
      SB=XY/R
      GO TO 50
   10 IF (Z2-Z1) 20,30,40
   20 CA=-1.D0
      CB=-1.D0
      SA=0.D0
      SB=0.D0
      GO TO 50
   30 CA=0.D0
      CB=0.D0
      SA=0.D0
      SB=0.D0
      GO TO 50
   40 CA=1.D0
      CB=1.D0
      SA=0.D0
      SB=0.D0
   50 CONTINUE
      CO=0
      DO 60 I=1,75
   60 C(I)=0.D0
      IF (PQ1.GT.PQ2) GO TO 70
      PQ=PQ2
      GO TO 80
   70 PQ=PQ1
   80 CONTINUE
      C(37)=1.D0
      IF (PQ.LT.2) GO TO 90
      C(56)=CA*CB
      C(41)=CA*SB
      C(26)=-SA
      C(53)=-SB
      C(38)=CB
      C(23)=0.D0
      C(50)=SA*CB
      C(35)=SA*SB
      C(20)=CA
      IF (PQ.LT.3) GO TO 90
      C2A=2*CA*CA-1.D0
      C2B=2*CB*CB-1.D0
      S2A=2*SA*CA
      S2B=2*SB*CB
      C(75)=C2A*CB*CB+0.5D0*C2A*SB*SB
      C(60)=0.5D0*C2A*S2B
      C(45)=0.8660254037841D0*C2A*SB*SB
      C(30)=-S2A*SB
      C(15)=-S2A*CB
      C(72)=-0.5D0*CA*S2B
      C(57)=CA*C2B
      C(42)=0.8660254037841D0*CA*S2B
      C(27)=-SA*CB
      C(12)=SA*SB
      C(69)=0.5773502691894D0*SB*SB*1.5D0
      C(54)=-0.8660254037841D0*S2B
      C(39)=CB*CB-0.5D0*SB*SB
      C(66)=-0.5D0*SA*S2B
      C(51)=SA*C2B
      C(36)=0.8660254037841D0*SA*S2B
      C(21)=CA*CB
      C(6)=-CA*SB
      C(63)=S2A*CB*CB+0.5D0*S2A*SB*SB
      C(48)=0.5D0*S2A*S2B
      C(33)=0.8660254037841D0*S2A*SB*SB
      C(18)=C2A*SB
      C(3)=C2A*CB
   90 CONTINUE
      RETURN
      END
      SUBROUTINE DANG(A1,A2,B1,B2,RCOS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         SAVE                                                           GL0892
**********************************************************************
*
*    DANG  DETERMINES THE ANGLE BETWEEN THE POINTS (A1,A2), (0,0),
*          AND (B1,B2).  THE RESULT IS PUT IN RCOS.
*
**********************************************************************
      PI=2.0D0* ASIN(1.0D00)
      ZERO=1.0D-6
      IF( ABS(A1).LT.ZERO.AND. ABS(A2).LT.ZERO) GO TO 10
      IF( ABS(B1).LT.ZERO.AND. ABS(B2).LT.ZERO) GO TO 10
      ANORM=1.0D0/ SQRT(A1**2+A2**2)
      BNORM=1.0D0/ SQRT(B1**2+B2**2)
      A1=A1*ANORM
      A2=A2*ANORM
      B1=B1*BNORM
      B2=B2*BNORM
      SINTH=(A1*B2)-(A2*B1)
      COSTH=A1*B1+A2*B2
      IF(COSTH.GT.1.0D0) COSTH=1.0D0
      IF(COSTH.LT.-1.0D0) COSTH=-1.0D0
      RCOS= ACOS(COSTH)
      IF( ABS(RCOS).LT.4.0D-4) GO TO 10
      IF(SINTH.GT.0.D0) RCOS=4.0D0* ASIN(1.0D00)-RCOS
      RCOS=-RCOS
      RETURN
   10 RCOS=0.0D0
      RETURN
      END
      SUBROUTINE DELMOL(COORD,I,J,NI,NJ,IA,ID,JA,JD,IX,RIJ,TOMB,ISP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION COORD(3,25)
      COMMON /DERIVS/ DS(16),DG(22),DR(100),TDX(3),TDY(3),TDZ(3)
      COMMON /EXTRA/  G(22),TX(3),TY(3),TZ(3)
         SAVE                                                           GL0892
      IF(NI.GT.1.OR.NJ.GT.1) CALL ROTAT(COORD,I,J,IX,RIJ,TOMB,2)
      IB=MAX(IA,ID)
      JB=MAX(JA,JD)
      DO 10 K=IA,IB
         KK=K-IA
         DO 10 L=K,IB
            LL=L-IA
            DO 10 M=JA,JB
               MM=M-JA
               DO 10 N=M,JB
                  NN=N-JA
                  ISP=ISP+1
                  IF(NN.EQ.0)THEN
                     IF(LL.EQ.0) THEN
C   (SS/SS)
                        DR(ISP)=DG(1)
                     ELSEIF(KK.EQ.0) THEN
C   (SP/SS)
                        DR(ISP)=DG(2)*TX(LL)+G(2)*TDX(LL)
                     ELSE
C   (PP/SS)
                        DR(ISP)=DG(3)*TX(KK)*TX(LL)
     1       +G(3)*(TDX(KK)*TX(LL)+TX(KK)*TDX(LL))
     2       +DG(4)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))
     3       +G(4)*(TDY(KK)*TY(LL)+TY(KK)*TDY(LL)
     4             +TDZ(KK)*TZ(LL)+TZ(KK)*TDZ(LL))
                     ENDIF
                  ELSEIF(MM.EQ.0) THEN
                     IF(LL.EQ.0) THEN
C   (SS/SP)
                        DR(ISP)=DG(5)*TX(NN)+G(5)*TDX(NN)
                     ELSEIF(KK.EQ.0) THEN
C   (SP/SP)
                        DR(ISP)=DG(6)*TX(LL)*TX(NN)
     1       +G(6)*(TDX(LL)*TX(NN)+TX(LL)*TDX(NN))
     2       +DG(7)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))
     3       +G(7)*(TDY(LL)*TY(NN)+TY(LL)*TDY(NN)
     4             +TDZ(LL)*TZ(NN)+TZ(LL)*TDZ(NN))
                     ELSE
C   (PP/SP)
                        DR(ISP)=DG(8)*TX(KK)*TX(LL)*TX(NN)
     1       +G(8)*(TDX(KK)*TX(LL)*TX(NN)+TX(KK)*TDX(LL)*TX(NN)
     2             +TX(KK)*TX(LL)*TDX(NN))
     3       +DG(9)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TX(NN)
     4       +G(9)*((TDY(KK)*TY(LL)+TY(KK)*TDY(LL)
     5              +TDZ(KK)*TZ(LL)+TZ(KK)*TDZ(LL))*TX(NN)
     6             +(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TDX(NN))
     7       +DG(10)*(TX(KK)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))
     8               +TX(LL)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN)))
     9       +G(10)*(TDX(KK)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))
     1              +TDX(LL)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))
     2              +TX(KK)*(TDY(LL)*TY(NN)+TY(LL)*TDY(NN)
     3                      +TDZ(LL)*TZ(NN)+TZ(LL)*TDZ(NN))
     4              +TX(LL)*(TDY(KK)*TY(NN)+TY(KK)*TDY(NN)
     5                      +TDZ(KK)*TZ(NN)+TZ(KK)*TDZ(NN)))
                     ENDIF
                  ELSEIF(LL.EQ.0) THEN
C   (SS/PP)
                     DR(ISP)=DG(11)*TX(MM)*TX(NN)
     1       +G(11)*(TDX(MM)*TX(NN)+TX(MM)*TDX(NN))
     2       +DG(12)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))
     3       +G(12)*(TDY(MM)*TY(NN)+TY(MM)*TDY(NN)
     4              +TDZ(MM)*TZ(NN)+TZ(MM)*TDZ(NN))
                  ELSEIF(KK.EQ.0) THEN
C   (SP/PP)
                     DR(ISP)=DG(13)*TX(LL)*TX(MM)*TX(NN)
     1       +G(13)*(TDX(LL)*TX(MM)*TX(NN)+TX(LL)*TDX(MM)*TX(NN)
     2              +TX(LL)*TX(MM)*TDX(NN))
     3       +DG(14)*TX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))
     4       +G(14)*(TDX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))
     5              +TX(LL)*(TDY(MM)*TY(NN)+TY(MM)*TDY(NN)
     6                      +TDZ(MM)*TZ(NN)+TZ(MM)*TDZ(NN)))
     7       +DG(15)*(TY(LL)*(TY(MM)*TX(NN)+TY(NN)*TX(MM))
     8               +TZ(LL)*(TZ(MM)*TX(NN)+TZ(NN)*TX(MM)))
     9       +G(15)*(TDY(LL)*(TY(MM)*TX(NN)+TY(NN)*TX(MM))
     1              +TDZ(LL)*(TZ(MM)*TX(NN)+TZ(NN)*TX(MM))
     2              +TY(LL)*(TDY(MM)*TX(NN)+TY(MM)*TDX(NN)
     3                      +TDY(NN)*TX(MM)+TY(NN)*TDX(MM))
     4              +TZ(LL)*(TDZ(MM)*TX(NN)+TZ(MM)*TDX(NN)
     5                      +TDZ(NN)*TX(MM)+TZ(NN)*TDX(MM)))
                  ELSE
C   (PP/PP)
                     DR(ISP)=DG(16)*TX(KK)*TX(LL)*TX(MM)*TX(NN)
     1       +G(16)*(TDX(KK)*TX(LL)*TX(MM)*TX(NN)
     2              +TX(KK)*TDX(LL)*TX(MM)*TX(NN)
     3              +TX(KK)*TX(LL)*TDX(MM)*TX(NN)
     4              +TX(KK)*TX(LL)*TX(MM)*TDX(NN))
     5       +DG(17)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TX(MM)*TX(NN)
     6       +G(17)*((TDY(KK)*TY(LL)+TY(KK)*TDY(LL)
     7               +TDZ(KK)*TZ(LL)+TZ(KK)*TDZ(LL))*TX(MM)*TX(NN)
     8              +(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))
     9              *(TDX(MM)*TX(NN)+TX(MM)*TDX(NN)))
     1       +DG(18)*TX(KK)*TX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))
     2       +G(18)*((TDX(KK)*TX(LL)+TX(KK)*TDX(LL))
     3                 *(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))
     4              +TX(KK)*TX(LL)*(TDY(MM)*TY(NN)+TY(MM)*TDY(NN)
     5                             +TDZ(MM)*TZ(NN)+TZ(MM)*TDZ(NN)))
                     DR(ISP)=DR(ISP)
     1       +DG(19)*(TY(KK)*TY(LL)*TY(MM)*TY(NN)
     2                  +TZ(KK)*TZ(LL)*TZ(MM)*TZ(NN))
     3       +G(19)*(TDY(KK)*TY(LL)*TY(MM)*TY(NN)
     4                 +TY(KK)*TDY(LL)*TY(MM)*TY(NN)
     5                 +TY(KK)*TY(LL)*TDY(MM)*TY(NN)
     6                 +TY(KK)*TY(LL)*TY(MM)*TDY(NN)
     7                 +TDZ(KK)*TZ(LL)*TZ(MM)*TZ(NN)
     8                 +TZ(KK)*TDZ(LL)*TZ(MM)*TZ(NN)
     9                 +TZ(KK)*TZ(LL)*TDZ(MM)*TZ(NN)
     1                 +TZ(KK)*TZ(LL)*TZ(MM)*TDZ(NN))
     2       +DG(20)*(TX(KK)*(TX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))
     3                          +TX(NN)*(TY(LL)*TY(MM)+TZ(LL)*TZ(MM)))
     4                  +TX(LL)*(TX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))
     5                          +TX(NN)*(TY(KK)*TY(MM)+TZ(KK)*TZ(MM))))
C      TO AVOID COMPILER DIFFICULTIES THIS IS DIVIDED
                     TEMP1=         TDX(KK)*(TX(MM)*(TY(LL)*TY(NN)+TZ(LL
     1)*TZ(NN))                          +TX(NN)*(TY(LL)*TY(MM)+TZ(LL)*T
     2Z(MM)))                 +TDX(LL)*(TX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(
     3NN))                          +TX(NN)*(TY(KK)*TY(MM)+TZ(KK)*TZ(MM)
     4))                 +TX(KK)*(TDX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))
     5                         +TDX(NN)*(TY(LL)*TY(MM)+TZ(LL)*TZ(MM)))
     6                 +TX(LL)*(TDX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))
     7                         +TDX(NN)*(TY(KK)*TY(MM)+TZ(KK)*TZ(MM)))
                     TEMP2=            TX(KK)*(TX(MM)*(TDY(LL)*TY(NN)+TY
     1(LL)*TDY(NN)                                 +TDZ(LL)*TZ(NN)+TZ(LL
     2)*TDZ(NN))                         +TX(NN)*(TDY(LL)*TY(MM)+TY(LL)*
     3TDY(MM)                                 +TDZ(LL)*TZ(MM)+TZ(LL)*TDZ
     4(MM)))                 +TX(LL)*(TX(MM)*(TDY(KK)*TY(NN)+TY(KK)*TDY(
     5NN)                                 +TDZ(KK)*TZ(NN)+TZ(KK)*TDZ(NN)
     6)                         +TX(NN)*(TDY(KK)*TY(MM)+TY(KK)*TDY(MM)
     7                                 +TDZ(KK)*TZ(MM)+TZ(KK)*TDZ(MM)))
                     DR(ISP)=DR(ISP)+G(20)*(TEMP1+TEMP2)
                     DR(ISP)=DR(ISP)
     1       +DG(21)*(TY(KK)*TY(LL)*TZ(MM)*TZ(NN)
     2                 +TZ(KK)*TZ(LL)*TY(MM)*TY(NN))
     3       +G(21)*(TDY(KK)*TY(LL)*TZ(MM)*TZ(NN)
     4                 +TY(KK)*TDY(LL)*TZ(MM)*TZ(NN)
     5                 +TY(KK)*TY(LL)*TDZ(MM)*TZ(NN)
     6                 +TY(KK)*TY(LL)*TZ(MM)*TDZ(NN)
     7                 +TDZ(KK)*TZ(LL)*TY(MM)*TY(NN)
     8                 +TZ(KK)*TDZ(LL)*TY(MM)*TY(NN)
     9                 +TZ(KK)*TZ(LL)*TDY(MM)*TY(NN)
     1                 +TZ(KK)*TZ(LL)*TY(MM)*TDY(NN))
                     DR(ISP)=DR(ISP)
     1       +DG(22)*(TY(KK)*TZ(LL)+TZ(KK)*TY(LL))
     2                 *(TY(MM)*TZ(NN)+TZ(MM)*TY(NN))
     3       +G(22)*((TDY(KK)*TZ(LL)+TY(KK)*TDZ(LL)
     4                  +TDZ(KK)*TY(LL)+TZ(KK)*TDY(LL))
     5                 *(TY(MM)*TZ(NN)+TZ(MM)*TY(NN))
     6                 +(TY(KK)*TZ(LL)+TZ(KK)*TY(LL))
     7                 *(TDY(MM)*TZ(NN)+TY(MM)*TDZ(NN)
     8                  +TDZ(MM)*TY(NN)+TZ(MM)*TDY(NN)))
                  ENDIF
   10 CONTINUE
      RETURN
      END
      SUBROUTINE DENROT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      COMMON /DENSTY/ P(MPACK),PA(MPACK),PB(MPACK)
      COMMON /GEOM  / GEO(3,NUMATM)
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /ELEMTS/ ELEMNT(107)
      COMMON /SCRACH/ B(MAXORB*MAXORB), BONDAB(MAXPAR**2-MAXORB*MAXORB)
************************************************************************
*
* DENROT PRINTS THE DENSITY MATRIX AS (S-SIGMA, P-SIGMA, P-PI) RATHER
*        THAN (S, PX, PY, PZ).
*
************************************************************************
      DIMENSION AROT(9,9), C(3,5,5), PAB(9,9), VECT(9,9)
      DIMENSION NATOM(MAXORB)
      DIMENSION XYZ(3,NUMATM), IROT(5,35), ISP(9)
      CHARACTER * 6 LINE(21)
      CHARACTER ELEMNT*2,ATORBS(9)*7,ITEXT(MAXORB)*7,JTEXT(MAXORB)*2
         SAVE                                                           GL0892
      DATA ATORBS/'S-SIGMA','P-SIGMA','  P-PI ','  P-PI ','D-SIGMA',
     1            '  D-PI ','  D-PI ',' D-DELL',' D-DELL'/
***********************************************************************
* IROT IS A MAPPING LIST. FOR EACH ELEMENT OF AROT 5 NUMBERS ARE
* NEEDED. THESE ARE, IN ORDER, FIRST AND SECOND SUBSCRIPTS OF AROT,
* AND FIRST,SECOND, AND THIRD SUBSCRIPTS OF C, THUS THE FIRST
* LINE OF IROT DEFINES AROT(1,1)=C(1,3,3)
*
***********************************************************************
      DATA IROT/1,1,1,3,3, 2,2,2,4,3, 3,2,2,2,3, 4,2,2,3,3, 2,3,2,4,2,
     1          3,3,2,2,2, 4,3,2,3,2, 2,4,2,4,4, 3,4,2,2,4, 4,4,2,3,4,
     2          5,5,3,1,5, 6,5,3,4,3, 7,5,3,3,3, 8,5,3,2,3, 9,5,3,5,3,
     3          5,6,3,1,2, 6,6,3,4,2, 7,7,3,3,2, 8,6,3,2,2, 9,6,3,5,2,
     4          5,7,3,1,4, 6,7,3,4,4, 7,7,3,3,4, 8,7,3,2,4, 9,7,3,5,4,
     5          5,8,3,1,1, 6,8,3,4,1, 7,8,3,3,1, 8,8,3,2,1, 9,8,3,5,1,
     6          5,9,3,1,5, 6,9,3,4,5, 7,9,3,3,5, 8,9,3,2,5, 9,9,3,5,5/
      DATA ISP /1,2,3,3,4,5,5,6,6/
      CALL GMETRY(GEO,XYZ)
      IPRT=0
      DO 120 I=1,NUMAT
         IF=NFIRST(I)
         IL=NLAST(I)
         IPQ=IL-IF-1
         II=IPQ+2
         IF(II.EQ.0)GOTO 120
         DO 10 I1=1,II
            J1=IPRT+ISP(I1)
            ITEXT(J1)=ATORBS(I1)
            JTEXT(J1)=ELEMNT(NAT(I))
            NATOM(J1)=I
   10    CONTINUE
         IPRT=J1
         IF(IPQ.NE.2)IPQ=MIN(MAX(IPQ,1),3)
         DO 110 J=1,I
            JF=NFIRST(J)
            JL=NLAST(J)
            JPQ=JL-JF-1
            JJ=JPQ+2
            IF(JJ.EQ.0)GOTO 110
            IF(JPQ.NE.2)JPQ=MIN(MAX(JPQ,1),3)
            DO 20 I1=1,9
               DO 20 J1=1,9
   20       PAB(I1,J1)=0.D0
            KK=0
            DO 30 K=IF,IL
               KK=KK+1
               LL=0
               DO 30 L=JF,JL
                  LL=LL+1
   30       PAB(KK,LL)=P(L+(K*(K-1))/2)
            CALL COE(XYZ(1,I),XYZ(2,I),XYZ(3,I),
     1                 XYZ(1,J),XYZ(2,J),XYZ(3,J),IPQ,JPQ,C,R)
            DO 40 I1=1,9
               DO 40 J1=1,9
   40       AROT(I1,J1)=0.D0
            DO 50 I1=1,35
   50       AROT(IROT(1,I1),IROT(2,I1))=
     1            C(IROT(3,I1),IROT(4,I1),IROT(5,I1))
            L1=ISP(II)
            L2=ISP(JJ)
            DO 60 I1=1,9
               DO 60 J1=1,9
   60       VECT(I1,J1)=-1.D0
            DO 70 I1=1,L1
               DO 70 J1=1,L2
   70       VECT(I1,J1)=0.D0
            IF(I.NE.J) THEN
               IJ=MAX(II,JJ)
               DO 90 I1=1,II
                  DO 90 J1=1,JJ
                     SUM=0.D0
                     DO 80 L1=1,IJ
                        DO 80 L2=1,IJ
   80                SUM=SUM+AROT(L1,I1)*PAB(L1,L2)*AROT(L2,J1)
   90          VECT(ISP(I1),ISP(J1))=
     1                        VECT(ISP(I1),ISP(J1))+SUM**2
            ENDIF
            K=0
            DO  100 I1=IF,IL
               K=K+1
               L=0
               DO 100 J1=JF,JL
                  L=L+1
  100       IF(J1.LE.I1) B(J1+(I1*(I1-1))/2)=VECT(K,L)
  110    CONTINUE
  120 CONTINUE
C
C NOW TO REMOVE ALL THE DEAD SPACE IN P, CHARACTERIZED BY -1.0
C
      LINEAR=(NORBS*(NORBS+1))/2
      L=0
      DO 130 I=1,LINEAR
         IF(B(I).GT.-0.1) THEN
            L=L+1
            B(L)=B(I)
         ENDIF
  130 CONTINUE
C
C   PUT ATOMIC ORBITAL VALENCIES ONTO THE DIAGONAL
C
      DO 160 I=1,IPRT
         SUM=0.D0
         II=(I*(I-1))/2
         DO 140 J=1,I
  140    SUM=SUM+B(J+II)
         DO 150 J=I+1,IPRT
  150    SUM=SUM+B((J*(J-1))/2+I)
  160 B((I*(I+1))/2)=SUM
      DO 170 I=1,21
  170 LINE(I)='------'
      LIMIT=(IPRT*(IPRT+1))/2
      KK=8
      NA=1
  180 LL=0
      M=MIN0((IPRT+1-NA),6)
      MA=2*M+1
      M=NA+M-1
      WRITE(6,'(/16X,10(1X,A7,3X))')(ITEXT(I),I=NA,M)
      WRITE(6,'(15X,10(2X,A2,I3,4X))')(JTEXT(I),NATOM(I),I=NA,M)
      WRITE (6,'(20A6)') (LINE(K),K=1,MA)
      DO 200 I=NA,IPRT
         LL=LL+1
         K=(I*(I-1))/2
         L=MIN0((K+M),(K+I))
         K=K+NA
         IF ((KK+LL).LE.50) GO TO 190
         WRITE (6,'(''1'')')
         WRITE(6,'(/17X,10(1X,A7,3X))')(ITEXT(N),N=NA,M)
         WRITE(6,'( 17X,10(2X,A2,I3,4X))')(JTEXT(N),NATOM(N),N=NA,M)
         WRITE (6,'(20A6)') (LINE(N),N=1,MA)
         KK=4
         LL=0
  190    WRITE (6,'(1X,A7,1X,A2,I3,10F11.6)')
     1   ITEXT(I),JTEXT(I),NATOM(I),(B(N),N=K,L)
  200 CONTINUE
      IF (L.GE.LIMIT) GO TO 210
      KK=KK+LL+4
      NA=M+1
      IF ((KK+IPRT+1-NA).LE.50) GO TO 180
      KK=4
      WRITE (6,'(''1'')')
      GO TO 180
  210 RETURN
      END
      SUBROUTINE DENSIT( C,MDIM, NORBS,NDUBL, NSINGL, FRACT, P, MODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION P(*), C(MDIM,*)
         SAVE                                                           GL0892
C***********************************************************************
C
C   DENSIT COMPUTES THE DENSITY MATRIX GIVEN THE EIGENVECTOR MATRIX, AND
C          INFORMATION ABOUT THE M.O. OCCUPANCY.
C
C  INPUT:  C     = SQUARE EIGENVECTOR MATRIX, C IS OF SIZE MDIM BY MDIM
C                  AND THE EIGENVECTORS ARE STORED IN THE TOP LEFT-HAND
C                  CORNER.
C          NORBS = NUMBER OF ORBITALS
C          NDUBL = NUMBER OF DOUBLY-OCCUPIED M.O.S ( =0 IF UHF)
C          NSINGL= NUMBER OF SINGLY OR FRACTIONALLY OCCUPIED M.O.S.
C
C   ON EXIT: P   = DENSITY MATRIX
C
C***********************************************************************
C
C SET UP LIMITS FOR SUMS
C  NL1 = BEGINING OF ONE ELECTRON SUM
C  NU1 = END OF SAME
C  NL2 = BEGINING OF TWO ELECTRON SUM
C  NU2 = END OF SAME
C
      NORBS2=NORBS/2
      IF(NSINGL .GT. NORBS2 .AND. MODE .EQ. 1) THEN
C
C    TAKE POSITRON EQUIVALENT
C
         SIGN=-1.D0
         FRAC=2.D0-FRACT
         IF(NDUBL.EQ.0)THEN
            CONST=1.D0
            NL2=2
            NU2=0
            NL1=NSINGL+1
            NU1=NORBS
         ELSE
            CONST=2.D0
            NL2=NSINGL+1
            NU2=NORBS
            NL1=NDUBL+1
            NU1=NSINGL
         ENDIF
      ELSE
C
C    TAKE ELECTRON EQUIVALENT
C
         SIGN=1.D0
         FRAC=FRACT
         CONST=0.D0
         NL2=1
         NU2=NDUBL
         NL1=NDUBL+1
         NU1=NSINGL
      ENDIF
      l = 1
      do 40 i=1,norbs
      l = l + i - 1
           do 30 j=l,l+i-1
           p(j) = 0.0
  30       continue
c
           do 20 k=nl2,nu2
                do 10 j=l,l+i-1
                p(j) = p(j) +(2.0*c(i,k))*c(j-l+1,k)
  10            continue
  20       continue
c
           do 50 k=nl1,nu1
                do 60 j=l,l+i-1
                p(j) = p(j) + (frac*c(i,k))*c(j-l+1,k)
  60            continue
  50       continue
c
           if(nsingl.gt.norbs2) then
                do 70 j=l,l+i-1
                p(j) = - p(j)
  70            continue
           p(l+i-1) = const + p(l+i-1)
           endif
c
  40  continue
      RETURN
      END
      SUBROUTINE DERIV(GEO,GRAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION GRAD(*), GEO(3,*)
      COMMON / EULER/ TVEC(3,3), ID
      COMMON /GEOVAR/ XDUMMY(MAXPAR), NVAR, LOC(2,MAXPAR)               IR0394 
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /GRAVEC/ COSINE
      COMMON /GEOSYM/ NDEP, IDUMYS(MAXPAR,3)
C    changed common path for portability  (IR)
      COMMON /PATHI / LATOM,LPARAM
      COMMON /PATHR / REACT(200)
c
      COMMON /UCELL / L1L,L2L,L3L,L1U,L2U,L3U
      COMMON /XYZGRA/ DXYZ(3,NUMATM*27)
      COMMON /ENUCLR/ ENUCLR
      COMMON /NUMCAL/ NUMCAL
      COMMON /DENSTY/ P(MPACK), PA(MPACK), PB(MPACK)
      COMMON /WMATRX/ WJ(N2ELEC), WK(N2ELEC)
      COMMON /HMATRX/ H(MPACK)
      COMMON /ATHEAT/ ATHEAT
C***********************************************************************
C
C    DERIV CALCULATES THE DERIVATIVES OF THE ENERGY WITH RESPECT TO THE
C          INTERNAL COORDINATES. THIS IS DONE BY FINITE DIFFERENCES.
C
C    THE MAIN ARRAYS IN DERIV ARE:
C        LOC    INTEGER ARRAY, LOC(1,I) CONTAINS THE ADDRESS OF THE ATOM
C               INTERNAL COORDINATE LOC(2,I) IS TO BE USED IN THE
C               DERIVATIVE CALCULATION.
C        GEO    ARRAY \GEO\ HOLDS THE INTERNAL COORDINATES.
C        GRAD   ON EXIT, CONTAINS THE DERIVATIVES
C
C***********************************************************************
      COMMON /KEYWRD / KEYWRD
      COMMON /ERRFN  / ERRFN(MAXPAR)
      CHARACTER*80 KEYWRD
      DIMENSION CHANGE(3), COORD(3,NUMATM), COLD(3,NUMATM*27)
     1,         XDERIV(3), XPARAM(MAXPAR), XJUC(3), W(N2ELEC)
      LOGICAL DEBUG, TIMES, HALFE, FAST, SCF1, CI, PRECIS
      EQUIVALENCE (W,WJ)
         SAVE                                                           GL0892
      DATA ICALCN /0/
      IF(ICALCN.NE.NUMCAL) THEN
         I = INDEX(KEYWRD,'PRESS')
         PRESS=0.D0
         IF(I.NE.0) PRESS=READA(KEYWRD,I)*1476.8992D0
         IDLO=NATOMS+1
         IF(LABELS(NATOMS) .EQ. 107) THEN
            IDLO=NATOMS
            IF(LABELS(NATOMS-1) .EQ. 107)THEN
               IDLO=NATOMS-1
               IF(LABELS(NATOMS-2) .EQ. 107)THEN
                  IDLO=NATOMS-2
               ENDIF
            ENDIF
         ENDIF
         GRLIM=0.01D0
         DEBUG = (INDEX(KEYWRD,'DERIV') .NE. 0)
         PRECIS= (INDEX(KEYWRD,'PRECIS') .NE. 0)
         TIMES = (INDEX(KEYWRD,'TIME') .NE. 0)
         CI    = (INDEX(KEYWRD,'C.I.') .NE. 0)
         SCF1  = (INDEX(KEYWRD,'1SCF') .NE. 0)
         ICALCN=NUMCAL
         IF(INDEX(KEYWRD,'RESTART') .EQ. 0) THEN
            DO 10 I=1,NVAR
   10       ERRFN(I)=0.D0
         ENDIF
         GRLIM=0.01D0
         IF(PRECIS)GRLIM=0.0001D0
         IF(INDEX(KEYWRD,'FULSCF') .GT.0) GRLIM=1.D9
         HALFE = (NOPEN.GT.NCLOSE .OR. CI)
         IDELTA=-7
*
*   IDELTA IS A MACHINE-PRECISION DEPENDANT INTEGER
*
         IF(HALFE.AND.PRECIS) IDELTA=-3
         IF(HALFE.AND..NOT.PRECIS) IDELTA=-3
         FAST=.TRUE.
         CHANGE(1)= 10.D0**IDELTA
         CHANGE(2)= 10.D0**IDELTA
         CHANGE(3)= 10.D0**IDELTA
C
C    CHANGE(I) IS THE STEP SIZE USED IN CALCULATING THE DERIVATIVES.
C    FOR "CARTESIAN" DERIVATIVES, CALCULATED USING DCART,AN
C    INFINITESIMAL STEP, HERE 0.000001, IS ACCEPTABLE. IN THE
C    HALF-ELECTRON METHOD A QUITE LARGE STEP IS NEEDED AS FULL SCF
C    CALCULATIONS ARE NEEDED, AND THE DIFFERENCE BETWEEN THE TOTAL
C    ENERGIES IS USED. THE STEP CANNOT BE VERY LARGE, AS THE SECOND
C    DERIVITIVE IN FLEPO IS CALCULATED FROM THE DIFFERENCES OF TWO
C    FIRST DERIVATIVES. CHANGE(1) IS FOR CHANGE IN BOND LENGTH,
C    (2) FOR ANGLE, AND (3) FOR DIHEDRAL.
C
         XDERIV(1)= 0.5D0/CHANGE(1)
         XDERIV(2)= 0.5D0/CHANGE(2)
         XDERIV(3)= 0.5D0/CHANGE(3)
      ENDIF
      GNORM=0.D0
      IF(NVAR.EQ.0) RETURN
      IF(DEBUG)THEN
         WRITE(6,'('' GEO AT START OF DERIV'')')
         WRITE(6,'(F19.5,2F12.5)')((GEO(J,I),J=1,3),I=1,NATOMS)
      ENDIF
      DO 20 I=1,NVAR
         XPARAM(I)=GEO(LOC(2,I),LOC(1,I))
   20 GNORM=GNORM+GRAD(I)**2
      GNORM=SQRT(GNORM)
      FAST=(GNORM .GT. GRLIM .AND. .NOT. SCF1 .OR. .NOT. HALFE)
cmgc      TIME1=SECOND()
      IF(NDEP.NE.0) CALL SYMTRY
      CALL GMETRY(GEO,COORD)
      IF( .NOT. FAST ) THEN
         IF(DEBUG)WRITE(6,'('' DOING FULL SCF''''S IN DERIV'')')
         CALL HCORE(COORD,H,W, WJ, WK, ENUCLR)
         IF(NORBS*NELECS.GT.0)THEN
            CALL ITER(H, W, WJ, WK, AA,.TRUE.,.FALSE.)
         ELSE
            AA=0.D0
         ENDIF
         LINEAR=(NORBS*(NORBS+1))/2
         DO 30 I=1,LINEAR
   30    P(I)=PA(I)*2.D0
         AA=(AA+ENUCLR)
      ENDIF
      CALL DCART(COORD,DXYZ)
      IF(NDEP.NE.0) CALL SYMTRY
      CALL GMETRY(GEO,COORD)
      IJ=0
      DO 70 II=1,NUMAT
         DO 60 IL=L1L,L1U
            DO 60 JL=L2L,L2U
               DO 60 KL=L3L,L3U
                  DO 40 LL=1,3
   40             XJUC(LL)=COORD(LL,II)+TVEC(LL,1)*IL+TVEC(LL,2)*JL+TVEC
     1(LL,3)*KL
                  IJ=IJ+1
                  DO 50 KK=1,3
                     COLD(KK,IJ)=XJUC(KK)
   50             CONTINUE
   60    CONTINUE
   70 CONTINUE
      SUM11=1.D-9
      SUM22=1.D-9
      SUM12=1.D-9
      DO 150 I=1,NVAR
         K=LOC(1,I)
         L=LOC(2,I)
         XSTORE=XPARAM(I)
         DO 80 J=1,NVAR
   80    GEO(LOC(2,J),LOC(1,J))=XPARAM(J)
         GEO(L,K)=XSTORE-CHANGE(L)
         IF(NDEP.NE.0) CALL SYMTRY
         CALL GMETRY(GEO,COORD)
C#         CALL GEOUT
C
C    USE LOOKUP TABLE OF CARTESIAN DERIVATIVES TO WORK OUT INTERNAL
C    COORDINATE DERIVATIVE.
C
         TOTL=0.D0
         IJ=0
         DO 130 II=1,NUMAT
            IF(ID.EQ.0) THEN
               DO 90 LL=1,3
   90          TOTL=TOTL+DXYZ(LL,II)*(COORD(LL,II)-COLD(LL,II))
            ELSE
               DO 120 IL=L1L,L1U
                  DO 120 JL=L2L,L2U
                     DO 120 KL=L3L,L3U
                        DO 100 LL=1,3
  100                   XJUC(LL)=COORD(LL,II)+TVEC(LL,1)*IL+TVEC(LL,2)*J
     1L+TVEC(LL,3)*KL
                        IJ=IJ+1
                        TOTL=TOTL+DXYZ(1,IJ)*(XJUC(1)-COLD(1,IJ))
     &                           +DXYZ(2,IJ)*(XJUC(2)-COLD(2,IJ))
     &                           +DXYZ(3,IJ)*(XJUC(3)-COLD(3,IJ))
  120          CONTINUE
            ENDIF
  130    CONTINUE
         TOTL=TOTL*XDERIV(L)
C
C   IF NEEDED, CALCULATE "EXACT" DERIVITIVES.
C
         IF( .NOT. FAST ) THEN
            CALL HCORE(COORD,H,W, WJ, WK,ENUCLR)
            IF(NORBS*NELECS.GT.0)THEN
               CALL ITER(H,W, WJ, WK,EE,.TRUE.,.FALSE.)
            ELSE
               EE=0.D0
            ENDIF
            DO 140 II=1,LINEAR
  140       P(II)=PA(II)*2.D0
            EE=(EE+ENUCLR)
            TOTL1=(AA-EE)*23.061D0*XDERIV(L)*2.D0
C#            WRITE(6,*)AA-EE
            ERRFN(I)=TOTL1-TOTL
         ENDIF
         GEO(L,K)=XSTORE
         SUM11=SUM11+GRAD(I)**2
         SUM22=SUM22+TOTL**2
         SUM12=SUM12+TOTL*GRAD(I)
         GRAD(I)=TOTL+ERRFN(I)
  150 CONTINUE
      IF(DEBUG) THEN
         WRITE(6,'('' GRADIENTS'')')
         WRITE(6,'(10F8.3)')(GRAD(I),I=1,NVAR)
         WRITE(6,'('' ERROR FUNCTION'')')
         WRITE(6,'(10F8.3)')(ERRFN(I),I=1,NVAR)
      ENDIF
      COSINE=SUM12/SQRT(SUM11*SUM22)
      IF(DEBUG)
     1WRITE(6,'('' COSINE OF SEARCH DIRECTION ='',F30.6)')COSINE
      IF( .NOT. FAST ) COSINE=1.D0
cmgc      IF(TIMES)
cmgc     1WRITE(6,'('' TIME FOR DERIVATIVES'',F12.6)')SECOND()-TIME1
      RETURN
      END
      SUBROUTINE DERS(M,N,RR,DEL1,DEL2,DEL3,IS,IOL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
************************************************************************
*                                                                      *
*    ON INPUT M    = INDEX OF FIRST ATOMIC ORBITAL                     *
*             N    = INDEX OF SECOND ATOMIC ORBITAL                    *
*             RR   = SQUARE IF INTERATOMIC DIATANCE (IN BOHR)          *
*             DEL1 = CATERSIAN DISTANCE IN DERIVATIVE DIRECTION        *
*             DEL2 = CARTESIAN DISTANCE IN M A.O.'S DIRECTION          *
*             DEL3 = CARTESIAN DISTANCE IN N A.O.'S DIRECTION          *
*             IS   = INDICATES TYPE OF A.O.-A.O. INTERACTION           *
*                  = 1 S/S, 2 S/P', 3 S/P, 4 P'/S, 5 P/S, 6 P/P',      *
*                    7 P'/P", 8 P'P', 9 P/P                            *
*             IOL  = INDEX FOR STORING DERIVATIVES IN DS               *
*                                                                      *
************************************************************************
      COMMON /DERIVS/ DS(16),DG(22),DR(100),TDX(3),TDY(3),TDZ(3)
      COMMON /TEMP/  CG(60,6),ZG(60,6)
      DIMENSION SS(6,6)
         SAVE                                                           GL0892
      A0=0.529167D0
      DO 110 I=1,6
         DO 110 J=1,6
            SS(I,J)=0.0D0
            APB=ZG(M,I)*ZG(N,J)
            AMB=ZG(M,I)+ZG(N,J)
            ADB=APB/AMB
            ADR=MIN(ADB*RR,35.D0)
            GO TO (10,20,30,40,50,60,70,80,90),IS
   10       ABN=-2.0D0*ADB*DEL1/(A0**2)
            GO TO 100
   20       ABN=-4.0D0*(ADB**2)*DEL1*DEL2/(SQRT(ZG(N,J))*(A0**3))
            GO TO 100
   30       ABN=(2.0D0*ADB/(SQRT(ZG(N,J))*A0))*
     1 (1.0D0-2.0D0*ADB*(DEL1**2)/(A0**2))
            GO TO 100
   40       ABN=4.0D0*(ADB**2)*DEL1*DEL2/(SQRT(ZG(M,I))*(A0**3))
            GO TO 100
   50       ABN=-(2.0D0*ADB/(SQRT(ZG(M,I))*A0))*
     1 (1.0D0-2.0D0*ADB*(DEL1**2)/(A0**2))
            GO TO 100
   60       ABN=-(4.0D0*(ADB**2)*DEL2/(SQRT(APB)*(A0**2)))*
     1 (1.0D0-2.0D0*ADB*(DEL1**2)/(A0**2))
            GO TO 100
   70       ABN=8.0D0*(ADB**3)*DEL1*DEL2*DEL3/(SQRT(APB)*(A0**4))
            GO TO 100
   80       ABN=-(8.0D0*(ADB**2)*DEL1/(SQRT(APB)*(A0**2)))*
     1 (0.5D0-ADB*(DEL2**2)/(A0**2))
            GO TO 100
   90       ABN=-(8.0D0*(ADB**2)*DEL1/(SQRT(APB)*(A0**2)))*
     1 (1.5D0-ADB*(DEL1**2)/(A0**2))
  100       SS(I,J)=SQRT((2.0D0*SQRT(APB)/AMB)**3)*EXP(-ADR)*ABN
  110 CONTINUE
      DO 120 I=1,6
         DO 120 J=1,6
            DS(IOL)=DS(IOL)+SS(I,J)*CG(M,I)*CG(N,J)
  120 CONTINUE
      RETURN
      END
      SUBROUTINE DHC (P,PA,PB,XI,NAT,IF,IM,IL,JF,JM,JL,
     1NORBS,DENER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION P(*), PA(*), PB(*)
      DIMENSION XI(3,*),NFIRST(2),NMIDLE(2),NLAST(2),NAT(*)
      dimension idiag(MAXORB),idiag1(MAXORB)
CGL0892  save idiag,idiag1
C***********************************************************************
C
C  DHC CALCULATES THE ENERGY CONTRIBUTIONS FROM THOSE PAIRS OF ATOMS
C         THAT HAVE BEEN MOVED BY SUBROUTINE DERIV.
C
C***********************************************************************
      COMMON /KEYWRD/ KEYWRD
     1       /ONELEC/ USS(107),UPP(107),UDD(107)
      COMMON /EULER / TVEC(3,3), ID
      COMMON /NUMCAL/ NUMCAL
      CHARACTER*80 KEYWRD
      LOGICAL UHF
      DIMENSION H(171), SHMAT(9,9), F(171),
     1          WJ(100), E1B(10), E2A(10), WK(100), W(100)
         SAVE                                                           GL0892
      DATA ICALCN /0/
      IF( ICALCN.NE.NUMCAL) THEN
         WLIM=4.D0
         IF(ID.EQ.0)WLIM=0.D0
         UHF=(INDEX(KEYWRD,'UHF') .NE. 0)
      ENDIF
      NFIRST(1)=1
      NMIDLE(1)=IM-IF+1
      NLAST(1)=IL-IF+1
      NFIRST(2)=NLAST(1)+1
      NMIDLE(2)=NFIRST(2)+JM-JF
      NLAST(2)=NFIRST(2)+JL-JF
      LINEAR=(NLAST(2)*(NLAST(2)+1))/2
      DO 10 I=1,LINEAR
         F(I)=0.D0
   10 H(I)=0.0D00
      if(ICALCN.NE.NUMCAL)  then
         do 11 i=1,maxorb
         idiag(i)=i*(i+1)/2
         idiag1(i)=i*(i-1)/2
11       continue
         ICALCN=NUMCAL
       endif
      JA=NFIRST(2)
      JB=NLAST(2)
      JC=NMIDLE(2)
      IA=NFIRST(1)
      IB=NLAST(1)
      IC=NMIDLE(1)
      JT=JB*(JB+1)/2
      J=2
      I=1
      NJ=NAT(2)
      NI=NAT(1)
      CALL H1ELEC(NI,NJ,XI(1,1),XI(1,2),SHMAT)
      IF(NAT(1).EQ.102.OR.NAT(2).EQ.102) THEN
         K=(JB*(JB+1))/2
         DO 30 J=1,K
   30    H(J)=0.D0
      ELSE
         J1=0
         DO 40 J=JA,JB
            jj=idiag1(j)
            J1=J1+1
            I1=0
            DO 40 I=IA,IB
               JJ=JJ+1
               I1=I1+1
               H(JJ)=SHMAT(I1,J1)
               F(JJ)=SHMAT(I1,J1)
   40    CONTINUE
      ENDIF
      KR=1
      IF(ID.EQ.0)THEN
         CALL ROTATE (NJ,NI,XI(1,2),XI(1,1),W(KR),KR,E2A,E1B,ENUCLR,100.
     1D0)
      ELSE
         CALL SOLROT (NJ,NI,XI(1,2),XI(1,1),WJ,WK,KR,E2A,E1B,ENUCLR,100.
     1D0)
      ENDIF
      IF(WJ(1).LT.WLIM)THEN
         DO 50 I=1,KR-1
   50    WK(I)=0.D0
      ENDIF
C
C    * ENUCLR IS SUMMED OVER CORE-CORE REPULSION INTEGRALS.
C
      I2=0
      DO 60 I1=IA,IC
         ii=idiag1(i1)+ia-1
         DO 60 J1=IA,I1
            II=II+1
            I2=I2+1
            H(II)=H(II)+E1B(I2)
   60 F(II)=F(II)+E1B(I2)
cdir$ ivdep
      DO  70 I1=IC+1,IB
         ii=idiag(i1)
         F(II)=F(II)+E1B(1)
   70 H(II)=H(II)+E1B(1)
      I2=0
      DO 80 I1=JA,JC
         ii=idiag1(i1)+ja-1
         DO 80 J1=JA,I1
            II=II+1
            I2=I2+1
            H(II)=H(II)+E2A(I2)
   80 F(II)=F(II)+E2A(I2)
cdir$ ivdep
      DO 90 I1=JC+1,JB
         ii=idiag(i1)
         F(II)=F(II)+E2A(1)
   90 H(II)=H(II)+E2A(1)
      CALL FOCK2D(F,P,PA,W, WJ, WK,2,NFIRST,NMIDLE,NLAST)
      EE=HELECT(NLAST(2),PA,H,F)
      IF( UHF ) THEN
         DO 100 I=1,LINEAR
  100    F(I)=H(I)
         CALL FOCK2D(F,P,PB,W, WJ, WK,2,NFIRST,NMIDLE,NLAST)
         EE=EE+HELECT(NLAST(2),PB,H,F)
      ELSE
         EE=EE*2.D0
      ENDIF
      DENER=EE+ENUCLR
      RETURN
C
      END
      FUNCTION DIAGI(IALPHA,IBETA,EIGA,XY,NMOS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION XY(NMECI,NMECI,NMECI,NMECI), EIGA(NMECI),
     1IALPHA(NMOS), IBETA(NMOS)
         SAVE                                                           GL0892
************************************************************************
*
*  CALCULATES THE ENERGY OF A MICROSTATE DEFINED BY IALPHA AND IBETA
*
************************************************************************
      X=0.0D0
      DO 20 I=1,NMOS
         IF (IALPHA(I).NE.0)THEN
            X=X+EIGA(I)
            DO 10  J=1,NMOS
               X=X+((XY(I,I,J,J)-XY(I,J,I,J))*IALPHA(J)*0.5D0 +
     1        (XY(I,I,J,J)            )*IBETA(J))
   10       CONTINUE
         ENDIF
   20 CONTINUE
      DO 40 I=1,NMOS
         IF (IBETA(I).NE.0) THEN
            X=X+EIGA(I)
            DO 30 J=1,I
   30       X=X+(XY(I,I,J,J)-XY(I,J,I,J))*IBETA(J)
         ENDIF
   40 CONTINUE
      DIAGI=X
      RETURN
      END
      SUBROUTINE DIAT(NI,NJ,XI,XJ,DI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
************************************************************************
*
*   DIAT CALCULATES THE DI-ATOMIC OVERLAP INTEGRALS BETWEEN ATOMS
*        CENTERED AT XI AND XJ.
*
*   ON INPUT NI  = ATOMIC NUMBER OF THE FIRST ATOM.
*            NJ  = ATOMIC NUMBER OF THE SECOND ATOM.
*            XI  = CARTESIAN COORDINATES OF THE FIRST ATOM.
*            XJ  = CARTESIAN COORDINATES OF THE SECOND ATOM.
*
*  ON OUTPUT DI  = DIATOMIC OVERLAP, IN A 9 * 9 MATRIX. LAYOUT OF
*                  ATOMIC ORBITALS IN DI IS
*                  1   2   3   4   5            6     7       8     9
*                  S   PX  PY  PZ  D(X**2-Y**2) D(XZ) D(Z**2) D(YZ)D(XY)
*
*   LIMITATIONS:  IN THIS FORMULATION, NI AND NJ MUST BE LESS THAN 107
*         EXPONENTS ARE ASSUMED TO BE PRESENT IN COMMON BLOCK EXPONT.
*
************************************************************************
      COMMON /KEYWRD/KEYWRD
      CHARACTER*80 KEYWRD
      INTEGER A,PQ2,B,PQ1,AA,BB
      LOGICAL ANALYT
      COMMON /EXPONT/ EMUS(107),EMUP(107),EMUD(107)
      DIMENSION DI(9,9),S(3,3,3),UL1(3),UL2(3),C(3,5,5),NPQ(107)
     1          ,XI(3),XJ(3), SLIN(27), IVAL(3,5)
     2, C1(3,5), C2(3,5), C3(3,5), C4(3,5), C5(3,5)
     3, S1(3,3), S2(3,3), S3(3,3)
      EQUIVALENCE(SLIN(1),S(1,1,1))
      EQUIVALENCE (C1(1,1),C(1,1,1)), (C2(1,1),C(1,1,2)),
     1            (C3(1,1),C(1,1,3)), (C4(1,1),C(1,1,4)),
     2            (C5(1,1),C(1,1,5)), (S1(1,1),S(1,1,1)),
     3            (S2(1,1),S(1,1,2)), (S3(1,1),S(1,1,3))
         SAVE                                                           GL0892
      DATA NPQ/1,0, 2,2,2,2,2,2,2,0, 0,3,3,3,3,3,3,0, 0,4,4,4,4,4,4,4,
     14,4,4,4,4,4,4,4,4,0, 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5
     2,32*6,15*0,3,5*0/
      DATA IVAL/1,0,9,1,3,8,1,4,7,1,2,6,0,0,5/
      ANALYT=(INDEX(KEYWRD,'ANALYT').NE.0)
      X1=XI(1)
      X2=XJ(1)
      Y1=XI(2)
      Y2=XJ(2)
      Z1=XI(3)
      Z2=XJ(3)
      PQ1=NPQ(NI)
      PQ2=NPQ(NJ)
      DO 20 I=1,9
         DO 10 J=1,9
            DI(I,J)=0.0D0
   10    CONTINUE
   20 CONTINUE
      CALL COE(X1,Y1,Z1,X2,Y2,Z2,PQ1,PQ2,C,R)
      IF(PQ1.EQ.0.OR.PQ2.EQ.0.OR.R.GE.10.D0) RETURN
      IF(R.LT.0.001)THEN
         RETURN
      ENDIF
      IA=MIN(PQ1,3)
      IB=MIN(PQ2,3)
      A=IA-1
      B=IB-1
      IF(ANALYT)THEN
         CALL GOVER(NI,NJ,XI,XJ,R,DI)
C#      WRITE(6,*)' OVERLAP FROM GOVER'
C#      WRITE(6,'(4F15.10)')SG
         RETURN
      ENDIF
      IF(NI.LT.18.AND.NJ.LT.18) THEN
         CALL DIAT2(NI,EMUS(NI),EMUP(NI),R,NJ,EMUS(NJ),EMUP(NJ),S)
      ELSE
         UL1(1)=EMUS(NI)
         UL2(1)=EMUS(NJ)
         UL1(2)=EMUP(NI)
         UL2(2)=EMUP(NJ)
         UL1(3)=MAX(EMUD(NI),0.3D0)
         UL2(3)=MAX(EMUD(NJ),0.3D0)
         DO 30 I=1,27
   30    SLIN(I)=0.0D0
         NEWK=MIN(A,B)
         NK1=NEWK+1
         DO 40 I=1,IA
            ISS=I
            IB=B+1
            DO 40 J=1,IB
               JSS=J
               DO 40 K=1,NK1
                  IF(K.GT.I.OR.K.GT.J) GOTO 40
                  KSS=K
                  S(I,J,K)=SS(PQ1,PQ2,ISS,JSS,KSS,UL1(I),UL2(J),R)       IR0295
   40    CONTINUE
      ENDIF
      DO 50 I=1,IA
         KMIN=4-I
         KMAX=2+I
         DO 50 J=1,IB
            IF(J.EQ.2)THEN
               AA=-1
               BB=1
            ELSE
               AA=1
               IF(J.EQ.3) THEN
                  BB=-1
               ELSE
                  BB=1
               ENDIF
            ENDIF
            LMIN=4-J
            LMAX=2+J
            DO 50 K=KMIN,KMAX
               DO 50 L=LMIN,LMAX
                  II=IVAL(I,K)
                  JJ=IVAL(J,L)
                  DI(II,JJ)=S1(I,J)*C3(I,K)*C3(J,L)*AA+
     1(C4(I,K)*C4(J,L)+C2(I,K)*C2(J,L))*BB*S2(I,J)+(C5(I,K)*C5(J,L)
     2+C1(I,K)*C1(J,L))*S3(I,J)
   50 CONTINUE
C#      WRITE(6,*)' OVERLAP FROM DIAT2'
C#      DO 12 I=1,4
C#  12  WRITE(6,'(4F15.10)')(DI(J,I),J=1,4)
      RETURN
      END
      SUBROUTINE DIAT2(NA,ESA,EPA,R12,NB,ESB,EPB,S)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION S(3,3,3)
C***********************************************************************
C
C OVERLP CALCULATES OVERLAPS BETWEEN ATOMIC ORBITALS FOR PAIRS OF ATOMS
C        IT CAN HANDLE THE ORBITALS 1S, 2S, 3S, 2P, AND 3P.
C
C***********************************************************************
      COMMON /SETC/ A(7),B(7),SA,SB,FACTOR,ISP,IPS
      DIMENSION INMB(17),III(78)
         SAVE                                                           GL0892
      DATA INMB/1,0,2,2,3,4,5,6,7,0,8,8,8,9,10,11,12/
C     NUMBERING CORRESPONDS TO BOND TYPE MATRIX GIVEN ABOVE
C      THE CODE IS
C
C     III=1      FIRST - FIRST  ROW ELEMENTS
C        =2      FIRST - SECOND
C        =3      FIRST - THIRD
C        =4      SECOND - SECOND
C        =5      SECOND - THIRD
C        =6      THIRD - THIRD
      DATA III/1,2,4,   2,4,4,   2,4,4,4,   2,4,4,4,4,
     1 2,4,4,4,4,4,   2,4,4,4,4,4,4,   3,5,5,5,5,5,5,6,
     2 3,5,5,5,5,5,5,6,6,   3,5,5,5,5,5,5,6,6,6,   3,5,5,5,5,5,5,6,6,6,6
     3, 3,5,5,5,5,5,5,6,6,6,6,6/
C
C      ASSIGNS BOND NUMBER
C
      JMAX=MAX0(INMB(NA),INMB(NB))
      JMIN=MIN0(INMB(NA),INMB(NB))
      NBOND=(JMAX*(JMAX-1))/2+JMIN
      II=III(NBOND)
      DO 10 I=1,3
         DO 10 J=1,3
            DO 10 K=1,3
   10 S(I,J,K)=0.D0
      RAB=R12/0.529167D0
      GOTO (20,30,40,50,60,70), II
C
C     ------------------------------------------------------------------
C *** THE ORDERING OF THE ELEMENTS WITHIN S IS
C *** S(1,1,1)=(S(B)/S(A))
C *** S(1,2,1)=(P-SIGMA(B)/S(A))
C *** S(2,1,1)=(S(B)/P-SIGMA(A))
C *** S(2,2,1)=(P-SIGMA(B)/P-SIGMA(A))
C *** S(2,2,2)=(P-PI(B)/P-PI(A))
C     ------------------------------------------------------------------
C *** FIRST ROW - FIRST ROW OVERLAPS
C
   20 CALL SET (ESA,ESB,NA,NB,RAB,NBOND,II)
      S(1,1,1)=.25D00*SQRT((SA*SB*RAB*RAB)**3)*(A(3)*B(1)-B(3)*A(1))
      RETURN
C
C *** FIRST ROW - SECOND ROW OVERLAPS
C
   30 CALL SET (ESA,ESB,NA,NB,RAB,NBOND,II)
      W=SQRT((SA**3)*(SB**5))*(RAB**4)*0.125D00
      S(1,1,1) = SQRT(1.D00/3.D00)
      S(1,1,1)=W*S(1,1,1)*(A(4)*B(1)-B(4)*A(1)+A(3)*B(2)-B(3)*A(2))
      IF (NA.GT.1) CALL SET (EPA,ESB,NA,NB,RAB,NBOND,II)
      IF (NB.GT.1) CALL SET (ESA,EPB,NA,NB,RAB,NBOND,II)
      W=SQRT((SA**3)*(SB**5))*(RAB**4)*0.125D00
      S(ISP,IPS,1)=W*(A(3)*B(1)-B(3)*A(1)+A(4)*B(2)-B(4)*A(2))
      RETURN
C
C *** FIRST ROW - THIRD ROW OVERLAPS
C
   40 CALL SET (ESA,ESB,NA,NB,RAB,NBOND,II)
      W=SQRT((SA**3)*(SB**7)/7.5D00)*(RAB**5)*0.0625D00
      SROOT3 = SQRT(3.D00)
      S(1,1,1)=W*(A(5)*B(1)-B(5)*A(1)+
     12.D00*(A(4)*B(2)-B(4)*A(2)))/SROOT3
      IF (NA.GT.1) CALL SET (EPA,ESB,NA,NB,RAB,NBOND,II)
      IF (NB.GT.1) CALL SET (ESA,EPB,NA,NB,RAB,NBOND,II)
      W=SQRT((SA**3)*(SB**7)/7.5D00)*(RAB**5)*0.0625D00
      S(ISP,IPS,1)=W*(A(4)*(B(1)+B(3))-B(4)*(A(1)+A(3))+
     1B(2)*(A(3)+A(5))-A(2)*(B(3)+B(5)))
      RETURN
C
C *** SECOND ROW - SECOND ROW OVERLAPS
C
   50 CALL SET (ESA,ESB,NA,NB,RAB,NBOND,II)
      W=SQRT((SA*SB)**5)*(RAB**5)*0.0625D00
      RT3=1.D00/SQRT(3.D00)
      S(1,1,1)=W*(A(5)*B(1)+B(5)*A(1)-2.0D00*A(3)*B(3))/3.0D00
      CALL SET (ESA,EPB,NA,NB,RAB,NBOND,II)
      IF (NA.GT.NB) CALL SET (EPA,ESB,NA,NB,RAB,NBOND,II)
      W=SQRT((SA*SB)**5)*(RAB**5)*0.0625D00
      D=A(4)*(B(1)-B(3))-A(2)*(B(3)-B(5))
      E=B(4)*(A(1)-A(3))-B(2)*(A(3)-A(5))
      S(ISP,IPS,1)=W*RT3*(D+E)
      CALL SET (EPA,ESB,NA,NB,RAB,NBOND,II)
      IF (NA.GT.NB) CALL SET (ESA,EPB,NA,NB,RAB,NBOND,II)
      W=SQRT((SA*SB)**5)*(RAB**5)*0.0625D00
      D=A(4)*(B(1)-B(3))-A(2)*(B(3)-B(5))
      E=B(4)*(A(1)-A(3))-B(2)*(A(3)-A(5))
      S(IPS,ISP,1)=-W*RT3*(E-D)
      CALL SET (EPA,EPB,NA,NB,RAB,NBOND,II)
      W=SQRT((SA*SB)**5)*(RAB**5)*0.0625D00
      S(2,2,1)=-W*(B(3)*(A(5)+A(1))-A(3)*(B(5)+B(1)))
      HD = .5D00
      S(2,2,2)=HD*W*(A(5)*(B(1)-B(3))-B(5)*(A(1)-A(3))
     1-A(3)*B(1)+B(3)*A(1))
      RETURN
C
C *** SECOND ROW - THIRD ROW OVERLAPS
C
   60 CALL SET (ESA,ESB,NA,NB,RAB,NBOND,II)
      W=SQRT((SA**5)*(SB**7)/7.5D00)*(RAB**6)*0.03125D00
      RT3 = 1.D00 / SQRT(3.D00)
      TD = 2.D00
      S(1,1,1)=W*(A(6)*B(1)+A(5)*B(2)-TD*(A(4)*B(3)+
     1A(3)*B(4))+A(2)*B(5)+A(
     21)*B(6))/3.D00
      CALL SET (ESA,EPB,NA,NB,RAB,NBOND,II)
      IF (NA.GT.NB) CALL SET (EPA,ESB,NA,NB,RAB,NBOND,II)
      W=SQRT((SA**5)*(SB**7)/7.5D00)*(RAB**6)*0.03125D00
      TD = 2.D00
      S(ISP,IPS,1)=W*RT3*(A(6)*B(2)+A(5)*B(1)-TD*(A(4)*B(4)+A(3)*B(3))
     1+A(2)*B(6)+A(1)*B(5))
      CALL SET (EPA,ESB,NA,NB,RAB,NBOND,II)
      IF (NA.GT.NB) CALL SET (ESA,EPB,NA,NB,RAB,NBOND,II)
      W=SQRT((SA**5)*SB**7/7.5D00)*(RAB**6)*0.03125D00
      TD = 2.D00
      S(IPS,ISP,1)=-W*RT3*(A(5)*(TD*B(3)-B(1))-B(5)*(TD*A(3)-A(1))-A(2
     1)*(B(6)-TD*B(4))+B(2)*(A(6)-TD*A(4)))
      CALL SET (EPA,EPB,NA,NB,RAB,NBOND,II)
      W=SQRT((SA**5)*SB**7/7.5D00)*(RAB**6)*0.03125D00
      S(2,2,1)=-W*(B(4)*(A(1)+A(5))-A(4)*(B(1)+B(5))
     1+B(3)*(A(2)+A(6))-A(3)*(B(2)+B(6)))
      HD = .5D00
      S(2,2,2)=HD*W*(A(6)*(B(1)-B(3))-B(6)*(A(1)-
     1A(3))+A(5)*(B(2)-B(4))-B(5
     2)*(A(2)-A(4))-A(4)*B(1)+B(4)*A(1)-A(3)*B(2)+B(3)*A(2))
      RETURN
C
C *** THIRD ROW - THIRD ROW OVERLAPS
C
   70 CALL SET (ESA,ESB,NA,NB,RAB,NBOND,II)
      W=SQRT((SA*SB*RAB*RAB)**7)/480.D00
      RT3 = 1.D00 / SQRT(3.D00)
      S(1,1,1)=W*(A(7)*B(1)-3.D00*(A(5)*B(3)-A(3)*B(5))-A(1)*B(7))/3.D00
      CALL SET (ESA,EPB,NA,NB,RAB,NBOND,II)
      IF (NA.GT.NB) CALL SET (EPA,ESB,NA,NB,RAB,NBOND,II)
      W=SQRT((SA*SB*RAB*RAB)**7)/480.D00
      D=A(6)*(B(1)-B(3))-2.D00*A(4)*(B(3)-B(5))+A(2)*(B(5)-B(7))
      E=B(6)*(A(1)-A(3))-2.D00*B(4)*(A(3)-A(5))+B(2)*(A(5)-A(7))
      S(ISP,IPS,1)=W*RT3*(D-E)
      CALL SET (EPA,ESB,NA,NB,RAB,NBOND,II)
      IF (NA.GT.NB) CALL SET (ESA,EPB,NA,NB,RAB,NBOND,II)
      W=SQRT((SA*SB*RAB*RAB)**7)/480.D00
      D=A(6)*(B(1)-B(3))-2.D00*A(4)*(B(3)-B(5))+A(2)*(B(5)-B(7))
      E=B(6)*(A(1)-A(3))-2.D00*B(4)*(A(3)-A(5))+B(2)*(A(5)-A(7))
      S(IPS,ISP,1)=-W*RT3*(-D-E)
      CALL SET (EPA,EPB,NA,NB,RAB,NBOND,II)
      W=SQRT((SA*SB*RAB*RAB)**7)/480.D00
      TD = 2.D00
      S(2,2,1)=-W*(A(3)*(B(7)+TD*B(3))-A(5)*(B(1)+
     1TD*B(5))-B(5)*A(1)+A(7)*B(3))
      HD = .5D00
      S(2,2,2)=HD*W*(A(7)*(B(1)-B(3))+B(7)*(A(1)-
     1A(3))+A(5)*(B(5)-B(3)-B(1)
     2)+B(5)*(A(5)-A(3)-A(1))+2.D00*A(3)*B(3))
      RETURN
C
      END
C     ******************************************************************
      DOUBLE PRECISION FUNCTION DIGIT(STRING,ISTART)
C     FORTRAN FUNCTION TO CONVERT NUMERIC FIELD TO DOUBLE PRECISION
C     NUMBER.  THE STRING IS ASSUMED TO BE CLEAN (NO INVALID DIGIT
C     OR CHARACTER COMBINATIONS FROM ISTART TO THE FIRST NONSPACE,
C     NONDIGIT, NONSIGN, AND NONDECIMAL POINT CHARACTER).
C
      CHARACTER STRING*(*)
      DOUBLE PRECISION C1,C2,DECIML
      LOGICAL SIGN
         SAVE                                                           GL0892
C
C     DEFINE ASCII VALUES OF NUMERIC FIELD CHARACTERS
      I0=ICHAR('0')
      I9=ICHAR('9')
      INEG=ICHAR('-')
      IPOS=ICHAR('+')
      IDOT=ICHAR('.')
      ISPC=ICHAR(' ')
C
      C1=0.D0
      C2=0.D0
      SIGN=.TRUE.
      L=LEN(STRING)
C
C     DETERMINE THE CONTRIBUTION TO THE NUMBER GREATER THAN ONE
      IDIG=0
      DO 10 I=ISTART,L
         N=ICHAR(STRING(I:I))
         IF(N.GE.I0.AND.N.LE.I9)THEN
            IDIG=IDIG+1
            C1=C1*1.D1+N-I0
         ELSEIF(N.EQ.INEG.OR.N.EQ.IPOS.OR.N.EQ.ISPC)THEN
            IF(N.EQ.INEG)SIGN=.FALSE.
         ELSEIF(N.EQ.IDOT)THEN
            GOTO 20
         ELSE
            GOTO 40
         ENDIF
   10 CONTINUE
C
C     DETERMINE THE CONTRIBUTION TO THE NUMBER LESS THAN THAN ONE
   20 DECIML=1.D0
      DO 30 J=I+1,L
         N=ICHAR(STRING(J:J))
         IF(N.GE.I0.AND.N.LE.I9)THEN
            DECIML=DECIML/1.D1
            C2=C2+(N-I0)*DECIML
         ELSEIF(N.NE.ISPC)THEN
            GOTO 40
         ENDIF
   30 CONTINUE
C
C     PUT THE PIECES TOGETHER
   40 DIGIT=C1+C2
      IF(.NOT.SIGN)DIGIT=-DIGIT
      RETURN
      END
      SUBROUTINE DIHED(XYZ,I,J,K,L,ANGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XYZ(3,*)
         SAVE                                                           GL0892
*********************************************************************
*
*      DIHED CALCULATES THE DIHEDRAL ANGLE BETWEEN ATOMS I, J, K,
*            AND L.  THE CARTESIAN COORDINATES OF THESE ATOMS
*            ARE IN ARRAY XYZ.
*
*     DIHED IS A MODIFIED VERSION OF A SUBROUTINE OF THE SAME NAME
*           WHICH WAS WRITTEN BY DR. W. THEIL IN 1973.
*
*********************************************************************
      XI1=XYZ(1,I)-XYZ(1,K)
      XJ1=XYZ(1,J)-XYZ(1,K)
      XL1=XYZ(1,L)-XYZ(1,K)
      YI1=XYZ(2,I)-XYZ(2,K)
      YJ1=XYZ(2,J)-XYZ(2,K)
      YL1=XYZ(2,L)-XYZ(2,K)
      ZI1=XYZ(3,I)-XYZ(3,K)
      ZJ1=XYZ(3,J)-XYZ(3,K)
      ZL1=XYZ(3,L)-XYZ(3,K)
C      ROTATE AROUND Z AXIS TO PUT KJ ALONG Y AXIS
      DIST= SQRT(XJ1**2+YJ1**2+ZJ1**2)
      COSA=ZJ1/DIST
      IF(COSA.GT.1.0D0) COSA=1.0D0
      IF(COSA.LT.-1.0D0) COSA=-1.0D0
      DDD=1.0D0-COSA**2
      IF(DDD.LE.0.0) GO TO 10
      YXDIST=DIST* SQRT(DDD)
      IF(YXDIST.GT.1.0D-6) GO TO 20
   10 CONTINUE
      XI2=XI1
      XL2=XL1
      YI2=YI1
      YL2=YL1
      COSTH=COSA
      SINTH=0.D0
      GO TO 30
   20 COSPH=YJ1/YXDIST
      SINPH=XJ1/YXDIST
      XI2=XI1*COSPH-YI1*SINPH
      XJ2=XJ1*COSPH-YJ1*SINPH
      XL2=XL1*COSPH-YL1*SINPH
      YI2=XI1*SINPH+YI1*COSPH
      YJ2=XJ1*SINPH+YJ1*COSPH
      YL2=XL1*SINPH+YL1*COSPH
C      ROTATE KJ AROUND THE X AXIS SO KJ LIES ALONG THE Z AXIS
      COSTH=COSA
      SINTH=YJ2/DIST
   30 CONTINUE
      YI3=YI2*COSTH-ZI1*SINTH
      YL3=YL2*COSTH-ZL1*SINTH
      CALL DANG(XL2,YL3,XI2,YI3,ANGLE)
      IF (ANGLE .LT. 0.) ANGLE=4.0D0* ASIN(1.0D00)+ANGLE
      IF (ANGLE .GE. 6.2831853D0 ) ANGLE=0.D0
      RETURN
      END
      SUBROUTINE DIPIND (DIPVEC)
C...............................................................
C  MODIFICATION OF DIPOLE SUBROUTINE FOR USE IN THE CALUCLATION
C  OF THE INDUCED DIPOLES FOR POLARIZABILITIES.
C...............................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      COMMON /CORE  / CORE(107)
      COMMON /DENSTY/ P(MPACK),PA(MPACK),PB(MPACK)
      COMMON /GEOM  / GEO(3,NUMATM)
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /NUMCAL/ NUMCAL
      COMMON /KEYWRD/ KEYWRD
      COMMON /ISTOPE/ AMS(107)
      COMMON /MULTIP/ DD(107), QQ(107), AM(107), AD(107), AQ(107)
      DIMENSION Q(MAXORB),Q2(MAXORB),DIPVEC(3),CENTER(3),
     1          COORD(3,NUMATM)
      CHARACTER*80 KEYWRD
C
C***********************************************************************
C     DIPOLE CALCULATES DIPOLE MOMENTS
C
C  ON INPUT P     = DENSITY MATRIX
C           Q     = TOTAL ATOMIC CHARGES, (NUCLEAR + ELECTRONIC)
C           NUMAT = NUMBER OF ATOMS IN MOLECULE
C           NAT   = ATOMIC NUMBERS OF ATOMS
C           NFIRST= START OF ATOM ORBITAL COUNTERS
C           COORD = COORDINATES OF ATOMS
C
C  OUTPUT  DIPOLE = DIPOLE MOMENT
C***********************************************************************
C
C     IN THE ZDO APPROXIMATION, ONLY TWO TERMS ARE RETAINED IN THE
C     CALCULATION OF DIPOLE MOMENTS.
C     1. THE POINT CHARGE TERM (INDEPENDENT OF PARAMETERIZATION).
C     2. THE ONE-CENTER HYBRIDIZATION TERM, WHICH ARISES FROM MATRIX
C     ELEMENTS OF THE FORM <NS/R/NP>. THIS TERM IS A FUNCTION OF
C     THE SLATER EXPONENTS (ZS,ZP) AND IS THUS DEPENDENT ON PARAMETER-
C     IZATION. THE HYBRIDIZATION FACTORS (HYF(I)) USED IN THIS SUB-
C     ROUTINE ARE CALCULATED FROM THE FOLLOWING FORMULAE.
C     FOR SECOND ROW ELEMENTS <2S/R/2P>
C     HYF(I)= 469.56193322*(SQRT(((ZS(I)**5)*(ZP(I)**5)))/
C           ((ZS(I) + ZP(I))**6))
C     FOR THIRD ROW ELEMENTS <3S/R/3P>
C     HYF(I)=2629.107682607*(SQRT(((ZS(I)**7)*(ZP(I)**7)))/
C           ((ZS(I) + ZP(I))**8))
C     FOR FOURTH ROW ELEMENTS AND UP :
C     HYF(I)=2*(2.10716)*DD(I)
C     WHERE DD(I) IS THE CHARGE SEPARATION IN ATOMIC UNITS
C
C
C     REFERENCES:
C     J.A.POPLE & D.L.BEVERIDGE: APPROXIMATE M.O. THEORY
C     S.P.MCGLYNN, ET AL: APPLIED QUANTUM CHEMISTRY
C
      DIMENSION DIP(4,3)
      DIMENSION HYF(107,2)
      LOGICAL FORCE, CHARGD
         SAVE                                                           GL0892
      DATA HYF(1,1)     / 0.0D00           /
      DATA   HYF(1,2) /0.0D0     /
      DATA   HYF(5,2) /6.520587D0/
      DATA   HYF(6,2) /4.253676D0/
      DATA   HYF(7,2) /2.947501D0/
      DATA   HYF(8,2) /2.139793D0/
      DATA   HYF(9,2) /2.2210719D0/
      DATA   HYF(14,2)/6.663059D0/
      DATA   HYF(15,2)/5.657623D0/
      DATA   HYF(16,2)/6.345552D0/
      DATA   HYF(17,2)/2.522964D0/
      DATA ICALCN /0/
C
C  SETUP FOR DIPOLE CALCULATION
C
      CALL CHRGE (P,Q2)
      DO 10 I = 1,NUMAT
         L = NAT(I)
         Q(I) = CORE(L) - Q2(I)
   10 CONTINUE
      CALL GMETRY (GEO,COORD)
C
      IF (ICALCN.NE.NUMCAL) THEN                                         IR0295
         DO 20 I=2,107
   20    HYF(I,1)= 5.0832*DD(I)
         WTMOL=0.D0
         SUM=0.D0
         DO 30 I=1,NUMAT
            WTMOL=WTMOL+AMS(NAT(I))
   30    SUM=SUM+Q(I)
         CHARGD=(ABS(SUM).GT.0.5D0)
         ICALCN=NUMCAL                                                   IR0295
         FORCE=(INDEX(KEYWRD,'FORCE') +INDEX(KEYWRD,'IRC').NE. 0)
         ITYPE=1
         IF(INDEX(KEYWRD,'MINDO') .NE. 0)ITYPE=2
      ENDIF
      IF(CHARGD)THEN
C
C   NEED TO RESET ION'S POSITION SO THAT THE CENTER OF MASS IS AT THE
C   ORIGIN.
C
         DO 40 I=1,3
   40    CENTER(I)=0.D0
         DO 50 I=1,3
            DO 50 J=1,NUMAT
   50    CENTER(I)=CENTER(I)+AMS(NAT(J))*COORD(I,J)
         DO 60 I=1,3
   60    CENTER(I)=CENTER(I)/WTMOL
         DO 70 I=1,3
            DO 70 J=1,NUMAT
   70    COORD(I,J)=COORD(I,J)-CENTER(I)
      ENDIF
      DO 80 I=1,4
         DO 80 J=1,3
   80 DIP(I,J)=0.0D00
      DO 100 I=1,NUMAT
         NI=NAT(I)
         IA=NFIRST(I)
         L=NLAST(I)-IA
         DO 90 J=1,L
            K=((IA+J)*(IA+J-1))/2+IA
   90    DIP(J,2)=DIP(J,2)-HYF(NI,ITYPE)*P(K)
         DO 100 J=1,3
  100 DIP(J,1)=DIP(J,1)+4.803D00*Q(I)*COORD(J,I)
      DO 110 J=1,3
  110 DIP(J,3)=DIP(J,2)+DIP(J,1)
      DO 120 J=1,3
  120 DIP(4,J)=SQRT(DIP(1,J)**2+DIP(2,J)**2+DIP(3,J)**2)
      DIPVEC(1)= -DIP(1,3)
      DIPVEC(2)= -DIP(2,3)
      DIPVEC(3)= -DIP(3,3)
      DIPCOM = DIP(4,3)
C      WRITE (6,60) ((DIP(I,J),I=1,4),J=1,3)
C   60 FORMAT (3(4F10.3))
      RETURN
C
      END
      SUBROUTINE ENPART(UHF,H,ALPHA,BETA,P,Q,COORD)
C----------------------------------------------------------*
C
C     NEW SUB. ENPART,  MODIFIED BY TSUNEO HIRANO 1986/6/3/
C
C---------------------------------------------------------*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      PARAMETER (NATMS2 = (NUMATM*(NUMATM+1))/2)
      DIMENSION H(*), ALPHA(*), BETA(*), P(*), Q(*), COORD(3,*)
C--- DEFINED HERE, AND TO BE USED FOR ENPART-PRINT ONLY ---*
      DIMENSION EX(NATMS2,3), EMAT(NUMATM,NUMATM)
C--- END OF DIMENSION DEFINITION ----------------- BY TH --*
      LOGICAL UHF, MINDO3, AM1
      CHARACTER*80 KEYWRD
      CHARACTER*2 ENAME(NUMATM),ELEMNT
C***********************************************************************
C
C *** ENERGY PARTITIONING WITHIN THE UMINDO/3 AND UMNDO SCHEME
C     ROUTINE WRITTEN BY S.OLIVELLA, BARCELONA NOV. 1979.
C     EXTENDED TO AM1 AND PM3 BY JJPS.
C
C   ON INPUT UHF     = .TRUE. IF A U.H.F. CALCULATION.
C            H       = ONE-ELECTRON MATRIX.
C            ALPHA   = ALPHA ELECTRON DENSITY.
C            BETA    = BETA ELECTRON DENSITY.
C            P       = TOTAL ELECTRON DENSITY.
C            Q       = ATOM ELECTRON DENSITIES.
C
C    NOTHING IS CHANGED ON EXIT.
C
C@***        ENAME   = ELEMENT NAME OF ATOM, I.E. C, CL, SI ETC.
C@
C***********************************************************************
      COMMON /ONELEC/ USS(107), UPP(107), UDD(107)
      COMMON /CORE  / CORE(107)
      COMMON /IDEAS / FN1(107,10),FN2(107,10),FN3(107,10)
      COMMON /ALPHA3/ ALP3(153)
      COMMON /TWOEL3/ F03(107)
      COMMON /ALPHA / ALP(107)
      COMMON /TWOELE/ GSS(107),GSP(107),GPP(107),GP2(107),HSP(107)
     1                ,GSD(107),GPD(107),GDD(107)
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
C@
      COMMON/ELEMTS/ELEMNT(107)
C@
      COMMON /WMATRX/ W(N2ELEC), DUMMY(N2ELEC)
      COMMON /KEYWRD/ KEYWRD
      PARAMETER (MDUMY=MAXPAR*MAXPAR-NUMATM*3-NATMS2*4)
      COMMON /SCRACH/ EA(NUMATM,2),EAT(NUMATM), E(NATMS2,4),
     1XDUMY(MDUMY)
         SAVE                                                           GL0892
      MINDO3=(INDEX(KEYWRD,'MINDO').NE.0)
C@ --------------------------*
C@    ASSIGN ENAME(I)
      DO 10 I=1,NUMAT
         IA=NAT(I)
   10 ENAME(I)=ELEMNT(IA)
C@    CLEAR EMAT(I,I)
      DO 20 I=1,NUMAT
   20 EMAT(I,I)=0.0D0
C@ --------------------------*
C *** RECALCULATE THE DENSITY MATRICES IN THE UHF SCHEME
C
      LINEAR=NORBS*(NORBS+1)/2
      IF( .NOT. UHF) THEN
         DO 30 I=1,LINEAR
   30    BETA(I)=ALPHA(I)
      ENDIF
C
C *** ONE-CENTER ENERGIES
      K=0
      DO 50 I=1,NUMAT
         IA=NFIRST(I)
         IB=NLAST(I)
         NI=NAT(I)
         EA(I,1)=0.0D0
         DO 40 J=IA,IB
            K=K+J
            T=UPP(NI)
            IF(J.EQ.IA) T=USS(NI)
   40    EA(I,1)=EA(I,1)+P(K)*T
         ISS=(IA*(IA+1))/2
         EA(I,2)=0.5D0*GSS(NI)*P(ISS)*P(ISS)
     1  -0.5D0*GSS(NI)*(ALPHA(ISS)*ALPHA(ISS)+BETA(ISS)*BETA(ISS))
         IF(IA.EQ.IB) GO TO 50
         IA1=IA+1
         IA2=IA+2
         IXX=IA1*IA2/2
         IYY=IA2*IB/2
         IZZ=(IB*(IB+1))/2
         IXY=IA1+IA2*IA1/2
         IXZ=IA1+IB*IA2/2
         IYZ=IA2+IB*IA2/2
         ISX=IA+IA1*IA/2
         ISY=IA+IA2*IA1/2
         ISZ=IA+IB*IA2/2
         SS1=P(IXX)*P(IXX)+P(IYY)*P(IYY)+P(IZZ)*P(IZZ)
         SS2=P(ISS)*(P(IXX)+P(IYY)+P(IZZ))
         SS3=P(IXX)*P(IYY)+P(IXX)*P(IZZ)+P(IYY)*P(IZZ)
         SS4=P(ISX)*P(ISX)+P(ISY)*P(ISY)+P(ISZ)*P(ISZ)
         SS5=P(IXY)*P(IXY)+P(IXZ)*P(IXZ)+P(IYZ)*P(IYZ)
         TT1=ALPHA(IXX)*ALPHA(IXX)+ALPHA(IYY)*ALPHA(IYY)
     1+ALPHA(IZZ)*ALPHA(IZZ)+BETA(IXX)*BETA(IXX)
     2+BETA(IYY)*BETA(IYY)+BETA(IZZ)*BETA(IZZ)
         TT2=ALPHA(ISS)*(ALPHA(IXX)+ALPHA(IYY)+ALPHA(IZZ))
     1   +BETA(ISS)*(BETA(IXX)+BETA(IYY)+BETA(IZZ))
         TT3=ALPHA(IXX)*ALPHA(IYY)+ALPHA(IXX)*ALPHA(IZZ)
     1+ALPHA(IYY)*ALPHA(IZZ)+BETA(IXX)*BETA(IYY)
     2+BETA(IXX)*BETA(IZZ)+BETA(IYY)*BETA(IZZ)
         TT4=ALPHA(ISX)*ALPHA(ISX)+ALPHA(ISY)*ALPHA(ISY)
     1+ALPHA(ISZ)*ALPHA(ISZ)+BETA(ISX)*BETA(ISX)
     2+BETA(ISY)*BETA(ISY)+BETA(ISZ)*BETA(ISZ)
         TT5=ALPHA(IXY)*ALPHA(IXY)+ALPHA(IXZ)*ALPHA(IXZ)
     1+ALPHA(IYZ)*ALPHA(IYZ)+BETA(IXY)*BETA(IXY)
     2+BETA(IXZ)*BETA(IXZ)+BETA(IYZ)*BETA(IYZ)
         EA(I,2)=EA(I,2)+0.5D0*GPP(NI)*SS1+GSP(NI)*SS2
     1+GP2(NI)*SS3+HSP(NI)*SS4*2.0D0+0.5D0*(GPP(NI)-GP2(NI))*SS5*2.0D0
     2                -0.5D0*GPP(NI)*TT1-GSP(NI)*TT4-GP2(NI)*TT5-
     3        HSP(NI)*(TT2+TT4)-0.5D0*(GPP(NI)-GP2(NI))*(TT3+TT5)
   50 CONTINUE
      AM1=(INDEX(KEYWRD,'AM1')+INDEX(KEYWRD,'PM3').NE.0)
      IF(MINDO3) THEN
         WRITE(6,'(///,10X,''TOTAL ENERGY PARTITIONING IN MINDO/3'')')
      ELSEIF( INDEX(KEYWRD,'PM3').NE.0 ) THEN
         WRITE(6,'(///,10X,''TOTAL ENERGY PARTITIONING IN PM3'')')
      ELSEIF( INDEX(KEYWRD,'AM1').NE.0 ) THEN
         WRITE(6,'(///,10X,''TOTAL ENERGY PARTITIONING IN AM1'')')
      ELSE
         WRITE(6,'(///,10X,''TOTAL ENERGY PARTITIONING IN MNDO'')')
      ENDIF
      KL=0
   60 K=KL+1
      KL=KL+10
      KL=MIN0(KL,NUMAT)
      DO 70 I=K,KL
   70 EAT(I)=EA(I,1)+EA(I,2)
      IF(NUMAT.GT.KL) GO TO 60
   80 EAU=0.0D0
      EAE=0.0D0
      DO 90 I=1,NUMAT
         EAU=EAU+EA(I,1)
   90 EAE=EAE+EA(I,2)
      TONE=EAU+EAE
C *** TWO-CENTER ENERGIES
C     RESONANCE (E(N,1)) TERMS
      N=1
      DO 110 II=2,NUMAT
         E(N,1)=0.0D0
         IA=NFIRST(II)
         IB=NLAST(II)
         IMINUS=II-1
         ONEII=1.D0
         IF(NAT(II).EQ.102)ONEII=0.D0
         DO 100 JJ=1,IMINUS
            N=N+1
            JA=NFIRST(JJ)
            JB=NLAST(JJ)
            ONEJJ=1.D0
            IF(NAT(JJ).EQ.102)ONEJJ=0.D0
            E(N,1)=0.0D0
            DO 100 I=IA,IB
               KA=(I*(I-1))/2
               DO 100 K=JA,JB
                  IK=KA+K
  100    E(N,1)=E(N,1)+2.0D0*P(IK)*H(IK)*ONEII*ONEJJ
  110 N=N+1
C
C     THE CODE THAT FOLLOWS APPLIES ONLY TO MNDO
C
      IF(.NOT.MINDO3) THEN
C     CORE-CORE REPULSION (E(N,2)) AND CORE-ELEC. ATTRACTION (E(N,3)).
         N=1
         KK=0
         DO 190 II=2,NUMAT
            E(N,2)=0.0D0
            E(N,3)=0.0D0
            IA=NFIRST(II)
            IB=NLAST(II)
            NI=NAT(II)
            ISS=(IA*(IA+1))/2
            IMINUS=II-1
            DO 180 JJ=1,IMINUS
               N=N+1
               JA=NFIRST(JJ)
               JB=NLAST(JJ)
               NJ=NAT(JJ)
               JSS=(JA*(JA+1))/2
               KK=KK+1
               G=W(KK)
               R=SQRT((COORD(1,II)-COORD(1,JJ))**2+(COORD(2,II)-COORD(2,
     1JJ))**2+  (COORD(3,II)-COORD(3,JJ))**2)
               SCALE=1.0D0+EXP(-ALP(NI)*R)+EXP(-ALP(NJ)*R)
               NT=NI+NJ
               IF(NT.LT.8.OR.NT.GT.9) GO TO 120
               IF(NI.EQ.7.OR.NI.EQ.8) SCALE=SCALE+(R-1.0D0)*EXP(-ALP(NI)
     1*R)
               IF(NJ.EQ.7.OR.NJ.EQ.8) SCALE=SCALE+(R-1.0D0)*EXP(-ALP(NJ)
     1*R)
  120          E(N,2)=CORE(NI)*CORE(NJ)*G*SCALE
               IF( AM1 )THEN
                  SCALE=0.0D0
                  DO 130 IG=1,10
                     IF(ABS(FN1(NI,IG)).GT.0.D0)
     1SCALE=SCALE +CORE(NI)*CORE(NJ)/R*
     2FN1(NI,IG)*EXP(-FN2(NI,IG)*(R-FN3(NI,IG))**2)
                     IF(ABS(FN1(NJ,IG)).GT.0.D0)
     1SCALE=SCALE +CORE(NI)*CORE(NJ)/R*
     2FN1(NJ,IG)*EXP(-FN2(NJ,IG)*(R-FN3(NJ,IG))**2)
  130             CONTINUE
                  E(N,2)=E(N,2)+SCALE
               ENDIF
               E(N,3)=-(P(ISS)*CORE(NJ)+P(JSS)*CORE(NI))*G
               IF(NJ.LT.3) GO TO 150
               KINC=9
               JAP1=JA+1
               DO 140 K=JAP1,JB
                  KC=(K*(K-1))/2
                  DO 140 L=JA,K
                     KL=KC+L
                     BB=2.0D0
                     IF(K.EQ.L) BB=1.0D0
                     KK=KK+1
  140          E(N,3)=E(N,3)-P(KL)*CORE(NI)*BB*W(KK)
               GO TO 160
  150          KINC=0
  160          IF(NI.LT.3) GO TO 180
               IAP1=IA+1
               DO 170 I=IAP1,IB
                  KA=(I*(I-1))/2
                  DO 170 J=IA,I
                     IJ=KA+J
                     AA=2.0D0
                     IF(I.EQ.J) AA=1.0D0
                     KK=KK+1
                     E(N,3)=E(N,3)-P(IJ)*CORE(NJ)*AA*W(KK)
  170          KK=KK+KINC
  180       CONTINUE
  190    N=N+1
C     COULOMB (E(N,4)) AND EXCHANGE (EX(N)) TERMS
         N=1
         KK=0
         DO 210 II=2,NUMAT
            E(N,4)=0.0D0
            EX(N,1)=0.0D0
            IA=NFIRST(II)
            IB=NLAST(II)
            IMINUS=II-1
            DO 200 JJ=1,IMINUS
               JA=NFIRST(JJ)
               JB=NLAST(JJ)
               N=N+1
               E(N,4)=0.0D0
               EX(N,1)=0.0D0
               DO 200 I=IA,IB
                  KA=(I*(I-1))/2
                  DO 200 J=IA,I
                     KB=(J*(J-1))/2
                     IJ=KA+J
                     AA=2.0D0
                     IF(I.EQ.J) AA=1.0D0
                     PIJ=P(IJ)
                     DO 200 K=JA,JB
                        KC=(K*(K-1))/2
                        IK=KA+K
                        JK=KB+K
                        DO 200 L=JA,K
                           IL=KA+L
                           JL=KB+L
                           KL=KC+L
                           BB=2.0D0
                           IF(K.EQ.L) BB=1.0D0
                           KK=KK+1
                           G=W(KK)
                           E(N,4)=E(N,4)+AA*BB*G*PIJ*P(KL)
  200       EX(N,1) = EX(N,1)
     1    -0.5D0*AA*BB*G*(ALPHA(IK)*ALPHA(JL)+ALPHA(IL)*ALPHA(JK)+
     2    BETA(IK)*BETA(JL)+BETA(IL)*BETA(JK))
  210    N=N+1
      ELSE
         N=1
         DO 270 I=2,NUMAT
            E(N,2)=0.0D0
            E(N,3)=0.0D0
            E(N,4)=0.0D0
            EX(N,1)=0.0D0
            IA=NFIRST(I)
            IB=NLAST(I)
            NI=NAT(I)
            IMINUS=I-1
            DO 260 J=1,IMINUS
               N=N+1
               JA=NFIRST(J)
               JB=NLAST(J)
               NJ=NAT(J)
               RIJ=(COORD(1,I)-COORD(1,J))**2+(COORD(2,I)-COORD(2,J))**2
     1+  (COORD(3,I)-COORD(3,J))**2
               GIJ=14.399D0/SQRT(RIJ+(7.1995D0/F03(NI)+7.1995D0/F03(NJ))
     1**2)
               PAB2=0.0D0
               IJ=MAX(NI,NJ)
               NBOND=(IJ*(IJ-1))/2+NI+NJ-IJ
               RIJ=SQRT(RIJ)
               IF(NBOND.EQ.22 .OR. NBOND .EQ. 29) GO TO 220
               GO TO 230
  220          SCALE=ALP3(NBOND)*EXP(-RIJ)
               GO TO 240
  230          SCALE=EXP(-ALP3(NBOND)*RIJ)
  240          CONTINUE
               E(N,2)=CORE(NI)*CORE(NJ)*GIJ+
     1     ABS(CORE(NI)*CORE(NJ)*(14.399D0/RIJ-GIJ)*SCALE)
               E(N,3)=(-Q(I)*CORE(NJ)-Q(J)*CORE(NI))*GIJ
               E(N,4)=Q(I)*Q(J)*GIJ
               DO 250 K=IA,IB
                  KK=(K*(K-1))/2
                  DO 250 L=JA,JB
                     LK=KK+L
  250          PAB2=PAB2+ALPHA(LK)*ALPHA(LK)+BETA(LK)*BETA(LK)
  260       EX(N,1) = -PAB2*GIJ
  270    N=N+1
      ENDIF
      NUMAT1=(NUMAT*(NUMAT+1))/2
      DO 280 I=1,4
  280 E(NUMAT1,I)=0.0D0
      DO 290 I=1,3
  290 EX(NUMAT1,I)=0.0D0
C@ --------------------------*
C-----PRINT OUT ONE AND TWO CENTER ENERGIES
C
C     E(I,1):     RESONANCE ENERGY
C     E(I,2):     NUCLEAR-NUCLEAR REPULSION ENERGY
C     E(I,3):     ELECTRON-NUCLEAR ATTRACTION ENERGY
C     E(I,4):     ELECTRON-ELECTRON REPULSION ENERGY
C     EX(I,1):    EXCHANGE  ENERGY
C     EX(I,2):    EXCHANGE + RESONANCE ENERGY
      WRITE(6,'(//,''       ONE AND TWO CENTER ENERGIES (EV) '')')
C
      WRITE(6,'(/,''  [RESONANCE TERM] (EV)'')')
      CALL VECPRT(E,NUMAT)
C
      WRITE(6,'(/,''  [EXCHANGE TERM] (EV)'')')
      CALL VECPRT(EX,NUMAT)
C
      WRITE(6,'(/,''  [RESONANCE + EXCHANGE] (EV)'')')
      DO 300 N=1,NUMAT1
  300 EX(N,2) =E(N,1) + EX(N,1)
C
C   ADD IN MONOCENTRIC EXCHANGE AND COULOMBIC TERM
C
      DO 310 I=1,NUMAT
  310 EX((I*(I+1))/2,2)=EA(I,2)
      CALL VECPRT(EX(1,2),NUMAT)
C
      WRITE(6,'(/,''  [ELECTRON - ELECTRON REPULSION] (EV)'')')
      CALL VECPRT(E(1,4),NUMAT)
C
      WRITE(6,'(/,''  [ELECTRON-NUCLEAR ATTRACTION] (EV)'')')
      DO 320 I=1,NUMAT
  320 E((I*(I+1))/2,3)=EA(I,1)
      CALL VECPRT(E(1,3),NUMAT)
C
      WRITE(6,'(/,''  [NUCLEAR-NUCLEAR REPULSION] (EV)'')')
      CALL VECPRT(E(1,2),NUMAT)
C
      DO 330 N=1,NUMAT1
  330 EX(N,3) =E(N,4) + E(N,3) + E(N,2)
C     PRINT OUT OF TOTAL COULOMB TERM
      WRITE(6,'(/,''  [TOTAL COULOMB TERM (E-E, E-N, AND N-N)] (EV)'')')
      CALL VECPRT(EX(1,3),NUMAT)
C     PRINT OUT OF TWO-CENTER SUM(OFF-DIAGONAL) +
C                  ONE-CENTER SUM(DIAGONAL).
      WRITE(6,'(/,''  [TWO-CENTER SUM (OFF-DIAGONAL), AND  '',
     1''ONE-CENTER SUM (DIAGONAL)] (EV)'')')
      DO 340 N=1,NUMAT1
  340 EX(N,3)=EX(N,3)+EX(N,2)
      CALL VECPRT(EX(1,3),NUMAT)
C@    TOTAL SUM
      EABR=0.0D0
      EABX=0.0D0
      EABEE=0.0D0
      EABEN=0.0D0
      EABNN=0.0D0
      DO 350 I=1,NUMAT
  350 E((I*(I+1))/2,3)=0.D0
      DO 360 I=1,NUMAT1
         EABR=EABR+E(I,1)
         EABX=EABX+EX(I,1)
         EABEE=EABEE+E(I,4)
         EABEN=EABEN+E(I,3)
         EABNN=EABNN+E(I,2)
  360 CONTINUE
      EABRX=EABR+EABX
      EABE=EABEE+EABEN+EABNN
      TTWO=EABRX+EABE
      ET=TONE+TTWO
C@ ***************************************************************
      WRITE(6,370)
  370 FORMAT(///,'***  SUMMARY OF ENERGY PARTITION  ***')
      WRITE(6,380)
  380 FORMAT(1H ,'---------------------------------------')
      WRITE(6,'(''     ONE-CENTER TERMS'')')
      WRITE(6,390) EAU
  390 FORMAT(/,' ELECTRON-NUCLEAR  (ONE-ELECTRON) ',F17.4,' EV')
      WRITE(6,400) EAE
  400 FORMAT(' ELECTRON-ELECTRON (TWO-ELECTRON) ',F17.4,' EV')
      WRITE(6,410) TONE
  410 FORMAT(/,' TOTAL OF ONE-CENTER TERMS        ',18X,F15.4,' EV')
      WRITE(6,380)
      WRITE(6,'(''     TWO-CENTER TERMS'')')
      WRITE(6,420) EABR
  420 FORMAT(/,' RESONANCE ENERGY',8X,F15.4,' EV')
      WRITE(6,430) EABX
  430 FORMAT(' EXCHANGE ENERGY ',8X,F15.4,' EV')
      WRITE(6,440) EABRX
  440 FORMAT(/,' EXCHANGE + RESONANCE ENERGY:       ',F15.4,' EV')
      WRITE(6,450) EABEE
  450 FORMAT(/,' ELECTRON-ELECTRON REPULSION',F12.4,' EV')
      WRITE(6,460) EABEN
  460 FORMAT(  ' ELECTRON-NUCLEAR ATTRACTION',F12.4,' EV')
      WRITE(6,470) EABNN
  470 FORMAT(  ' NUCLEAR-NUCLEAR REPULSION  ',F12.4,' EV')
      WRITE(6,480) EABE
  480 FORMAT(/,' TOTAL ELECTROSTATIC INTERACTION    ',F15.4,' EV',/)
      WRITE(6,490) TTWO
  490 FORMAT(' GRAND TOTAL OF TWO-CENTER TERMS   ',17X,F15.4,' EV')
      WRITE(6,380)
      WRITE(6,500) ET
  500 FORMAT(' ETOT (EONE + ETWO)   ',30X,F15.4,' EV'//)
      RETURN
      END
      SUBROUTINE FFHPOL (HEAT0,ATPOL,DIPVEC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      CHARACTER*1 AXIS(3)
      CHARACTER*80 KEYWRD
      LOGICAL LARGE,POLDIP,DEBUG
C***********************************************************************
C  SUBROUTINE FOR THE FINITE FIELD CALCULATION OF ELECTRIC RESPONSE
C  PROPERTIES (DIPOLE MOMENT, POLARIZABILITY, AND 1ST AND 2ND
C  HYPERPOLARIZABILITY.
C
C  HENRY A. KURTZ, DEPARTMENT OF CHEMISTRY
C                  MEMPHIS STATE UNIVERSITY
C                  MEMPHIS, TN   38152
C
C***********************************************************************
      COMMON /CORE  / CORE(107)
      COMMON /GEOM  / GEO(3,NUMATM)
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /COORD / COORD(3,NUMATM)
      COMMON /KEYWRD/ KEYWRD
      COMMON /FIELD / EFIELD(3)
      COMMON /EULER / TVEC(3,3),IDTVEC
C
C
C     DIPE4 AND DIPDP HOLD THE CALCULATED DIPOLE MOMENTS
C
C     APOLE4 AND APOLDP HOLD THE POLARIZABILITY TENSOR AS
C                                A PACKED ARRAY XX,XY,YY,XZ,YZ,ZZ
C
C     BETAE4 AND BETAEP HOLD THE FIRST HYPERPOLARIZABILITY
C                                1. XXX
C                                2. YYY     6. YXX
C                                3. ZZZ     7. YZZ
C                                4. XYY     8. ZXX
C                                5. XZZ     9. ZYY
C
      DIMENSION HEATE(3,2),DIPVEC(3),EIGS(3),VECTRS(9),
     1          DIPE4(3),APOLE4(6),BETAE4(9),GAMME4(6),
     2          DIPDP(3),APOLDP(6),BETADP(9),GAMMDP(6),
     3          DIP1P(3),DIP1M(3),DIP2P(3),DIP2M(3)
      DIMENSION IPTBD(6)
         SAVE                                                           GL0892
C
C     Initialize HEAT3M and HEAT3P to eliminate warning message
C
      DATA HEAT3M,HEAT3P /0.0D+00,0.0D+00/
C
      DATA IPTBD /5,7,4,9,6,8/
C Energy: a.u. to kcal/mole
      AUTOKC = 23.061D+00*27.2107D+00
C Length: a.u. to Angstrom
      AUTOA  = 0.529177D+00
C Dipole: a.u. to debye
      AUTODB = 2.541563D+00
C Electric Field: a.u. to volt/meter
      AUTOVM = 51.4257D+00
      NBDIP = 1
      NGDIP = 4
      NBCNT = 4
      NGCNT = 4
C
      DATA AXIS/'X','Y','Z'/
      LARGE = (INDEX(KEYWRD,'LARGE').NE.0)
      DEBUG = (INDEX(KEYWRD,'DEBUG').NE.0)
C
C  FIELD STRENGTH IN A.U.
C
      EFVAL=0.001D0
      IDIP=1
C#      READ (7,10) EFVAL,IDIP
   10 FORMAT (F10.5,I5)
      WRITE (6,20) EFVAL
   20 FORMAT (//' APPLIED ELECTRIC FIELD MAGNITUDE: ',F15.5)
      POLDIP = .FALSE.
      IF (IDIP.NE.0) POLDIP = .TRUE.
      SFE = 1.D00/EFVAL
      WRITE (6,30) 6.74834*ATPOL
   30 FORMAT (//' ATOMIC CONTRIBUTION TO THE POLARIZABILITY: ',F15.6)
C.......................................................................
C  CALCULATE THE POLARIZABILITY AND HYPERPOLARIZABILITIES ALONG
C  THE THREE PRINCIPLE AXES.  (THESE AXES DEPEND ON YOUR ARBITRARY
C  ORIENTATION AND MAY NOT BE THE TRUE PRINCIPLE AXES.)
C.......................................................................
      DO 160 ID = 1,3
         IF (DEBUG) THEN
            WRITE (6,40) AXIS(ID)
   40       FORMAT (//,' ****** ',A1,' DIRECTION *****',/)
         ENDIF
C
C ZERO THE FIELD
C
         DO 50 I = 1,3
            EFIELD(I) = 0.0D00
   50    CONTINUE
         HNUC = 0.0D00
         DO 60 I = 1,NUMAT
            HNUC = HNUC + EFVAL*GEO(ID,I)*CORE(NAT(I))*AUTOVM
   60    CONTINUE
         HNUC = HNUC*23.061D00
C +E(ID)
         EFIELD(ID) = EFVAL
         CALL COMPFG(GEO,.TRUE.,HEAT1P,.TRUE.,GRAD,.FALSE.)
         CALL DIPIND (DIP1P)
         DIIP = DIP1P(ID)
C -E(ID)
         EFIELD(ID) = -EFVAL
         CALL COMPFG(GEO,.TRUE.,HEAT1M,.TRUE.,GRAD,.FALSE.)
         CALL DIPIND (DIP1M)
         DIIM = DIP1M(ID)
C +2E(ID)
         EFIELD(ID) = 2.0D00*EFVAL
         CALL COMPFG(GEO,.TRUE.,HEAT2P,.TRUE.,GRAD,.FALSE.)
         CALL DIPIND (DIP2P)
C -2E(ID)
         EFIELD(ID) = -2.0D00*EFVAL
         CALL COMPFG(GEO,.TRUE.,HEAT2M,.TRUE.,GRAD,.FALSE.)
         CALL DIPIND (DIP2M)
C
C  CORRECT FOR ELECTRIC FIELD - NUCLEAR INTERACTIONS
C
         HEAT1P = HEAT1P + HNUC
         HEATE(ID,1) = HEAT1P
         HEAT1M = HEAT1M - HNUC
         HEATE(ID,2) = HEAT1M
         HEAT2P = HEAT2P + HNUC*2.D00
         HEAT2M = HEAT2M - HNUC*2.D00
C
         IF (DEBUG) THEN
            WRITE (6,70)
   70       FORMAT (' ENERGIES AT: ',5X,'F',21X,'2F',21X,'3F',/)
            WRITE (6,80) HEAT1P,HEAT2P,HEAT3P,HEAT1M,HEAT2M,HEAT3M
   80       FORMAT ('   + ',3(F20.10,3X),/,'   - ',3(F20.10,3X))
         ENDIF
C
C DIPOLE
C
         ETERM = (1.0D00/12.D00)*(HEAT2P - HEAT2M)
     1        - (2.0D00/3.0D00)*(HEAT1P - HEAT1M)
         DIPE4(ID) = ETERM*SFE/AUTOKC
C
C ALPHA
C
         IVL = (ID*(ID+1))/2
         ETERM = 2.5D00*HEAT0 - (4.D00/3.D00)*(HEAT1P + HEAT1M)
     1        + (1.D00/12.0D00)*(HEAT2P + HEAT2M)
         APOLE4(IVL) = ETERM*SFE*SFE/AUTOKC + ATPOL*6.74834
C
C BETA
C
         ETERM = (HEAT1P - HEAT1M) - 0.5D00*(HEAT2P - HEAT2M)
         BETAE4(ID) = ETERM*SFE*SFE*SFE/AUTOKC
C
C GAMMA
C
         ETERM = 4.0D00*(HEAT1P + HEAT1M) - (HEAT2P + HEAT2M)
     1        - 6.0D00*HEAT0
         GAMME4(ID) = ETERM*SFE*SFE*SFE*SFE/AUTOKC
C
C DIPOLE CALCULATIONS
C
         DMU = (2.0D00/3.0D00)*(DIP1P(ID) + DIP1M(ID))
     1      - (1.D00/6.0D00)*(DIP2P(ID) + DIP2M(ID))
         DIPDP(ID) = DMU/AUTODB
         AE = (2.0D00/3.0D00)*(DIP1P(ID) - DIP1M(ID))
     1     - (1.0D00/12.D00)*(DIP2P(ID) - DIP2M(ID))
         APOLDP(IVL) = AE*SFE/AUTODB
         BE = (1.D00/3.0D00)*(DIP2P(ID) + DIP2M(ID)
     1                     - DIP1P(ID) - DIP1M(ID))
         BETADP(ID) = BE*SFE*SFE/AUTODB
         GE = 0.5D00*(DIP2P(ID) - DIP2M(ID))
     1     - (DIP1P(ID) - DIP1M(ID))
         GAMMDP(ID) = GE*SFE*SFE*SFE/AUTODB
         DO 90 KD = 1,3
            IF (KD.LT.ID) THEN
               KVL = (ID*(ID-1))/2 + KD
               AKI = (2.0D00/3.0D00)*(DIP1P(KD) - DIP1M(KD))
     1         - (1.0D00/12.0D00)*(DIP2P(KD) - DIP2M(KD))
               APOLDP(KVL) = AKI*SFE/AUTODB
            ENDIF
            IF (KD.NE.ID) THEN
               BKII = (1.0D00/3.0D00)*(DIP2P(KD) + DIP2M(KD)
     1                           - DIP1P(KD) - DIP1M(KD))
               NBD = IPTBD(NBDIP)
               BETADP(NBD) = BKII*SFE*SFE/AUTODB
               NBDIP = NBDIP + 1
            ENDIF
   90    CONTINUE
C.......................................................................
C
C  NOW CALCULATE THE OFF AXIS RESULTS.
C
C.......................................................................
         IDM1 = ID - 1
         DO 150 JD = 1,IDM1
            HNUCJ = 0.0D00
            DO 100 I = 1,NUMAT
               HNUCJ = HNUCJ + EFVAL*GEO(JD,I)*CORE(NAT(I))*51.4257
  100       CONTINUE
            HNUCJ = HNUCJ*23.061
            DO 110 I = 1,3
               EFIELD(I) = 0.0D00
  110       CONTINUE
C
C DIAGONAL FIELDS WITH COMPONENTS EQUAL TO EFVAL
C
            EFIELD(ID) = EFVAL
            EFIELD(JD) = EFVAL
            CALL COMPFG(GEO,.TRUE.,HPP,.TRUE.,GRAD,.FALSE.)
            CALL DIPIND (DIP1P)
            DPP = DIP1P(ID)
            EFIELD(JD) = -EFVAL
            CALL COMPFG(GEO,.TRUE.,HPM,.TRUE.,GRAD,.FALSE.)
            CALL DIPIND (DIP1P)
            DPM = DIP1P(ID)
            EFIELD(ID) = -EFVAL
            CALL COMPFG(GEO,.TRUE.,HMM,.TRUE.,GRAD,.FALSE.)
            CALL DIPIND (DIP1P)
            DMM = DIP1P(ID)
            EFIELD(JD) = EFVAL
            CALL COMPFG(GEO,.TRUE.,HMP,.TRUE.,GRAD,.FALSE.)
            CALL DIPIND (DIP1P)
            DMP = DIP1P(ID)
            HPP = HPP + HNUC + HNUCJ
            HPM = HPM + HNUC - HNUCJ
            HMM = HMM - HNUC - HNUCJ
            HMP = HMP - HNUC + HNUCJ
            IF (DEBUG) THEN
               WRITE (6,120)
  120          FORMAT (/,' ',12X,'+,+',15X,'+,-',15X,'-,+',15X,'-,-')
               WRITE (6,130) HPP,HPM,HMP,HMM
  130          FORMAT ('  E ',4F15.6)
            ENDIF
C
C  DIAGONAL FIELDS WITH COMPONENTS EQUAL TO 2*EFVAL
C
            EFIELD(ID) = EFVAL*2.D00
            EFIELD(JD) = EFVAL*2.D00
            CALL COMPFG(GEO,.TRUE.,H2PP,.TRUE.,GRAD,.FALSE.)
            EFIELD(JD) = -EFVAL*2.D00
            CALL COMPFG(GEO,.TRUE.,H2PM,.TRUE.,GRAD,.FALSE.)
            EFIELD(ID) = -EFVAL*2.D00
            CALL COMPFG(GEO,.TRUE.,H2MM,.TRUE.,GRAD,.FALSE.)
            EFIELD(JD) = EFVAL*2.D00
            CALL COMPFG(GEO,.TRUE.,H2MP,.TRUE.,GRAD,.FALSE.)
            H2PP = H2PP + 2.0D00*(HNUC + HNUCJ)
            H2PM = H2PM + 2.0D00*(HNUC - HNUCJ)
            H2MM = H2MM - 2.0D00*(HNUC + HNUCJ)
            H2MP = H2MP - 2.0D00*(HNUC - HNUCJ)
            IF (DEBUG) THEN
               WRITE (6,140) H2PP,H2PM,H2MP,H2MM
  140          FORMAT (' 2E ',4F15.6)
            ENDIF
C
            ATERM = (1.0D00/48.0D00)*(H2PP - H2PM - H2MP + H2MM)
     1          - (1.0D00/3.0D00)*(HPP - HPM - HMP + HMM)
            AIJ = ATERM*SFE*SFE/AUTOKC
            IVL = (ID*(ID-1))/2 + JD
            APOLE4(IVL) = AIJ
            BTERM = 0.5D00*(HMM - HPP + HPM - HMP)
     1          + HEATE(JD,1) - HEATE(JD,2)
            BJII = BTERM*SFE*SFE*SFE/AUTOKC
            BETAE4(NBCNT) = BJII
            NBCNT = NBCNT + 1
            BTERM = 0.5D00*(HMM - HPP + HMP - HPM)
     1          + HEATE(ID,1) - HEATE(ID,2)
            BIJJ = BTERM*SFE*SFE*SFE/AUTOKC
            BETAE4(NBCNT) = BIJJ
            NBCNT = NBCNT + 1
C
            GTERM = -(HPP + HMM + HPM + HMP) - 4.0D00*HEAT0
     1           + 2.0D00*(HEATE(ID,1) + HEATE(ID,2))
     2           + 2.0D00*(HEATE(JD,1) + HEATE(JD,2))
            GIIJJ = GTERM*SFE*SFE*SFE*SFE/AUTOKC
            GAMME4(NGCNT) = GIIJJ
            GDIP = 0.5D00*(DPP - DMP + DPM - DMM) - (DIIP - DIIM)
            GAMMDP(NGCNT) = GDIP*SFE*SFE*SFE/AUTODB
            NGCNT = NGCNT + 1
  150    CONTINUE
C
  160 CONTINUE
C-----------------------------------------------------------------------
C  SUMMARIZE THE RESULTS
C-----------------------------------------------------------------------
      WRITE (6,170)
  170 FORMAT (//,' ',30('*'),' DIPOLE ',30('*'),//)
      DIPE4T = SQRT(DIPE4(1)*DIPE4(1) + DIPE4(2)*DIPE4(2)
     1              + DIPE4(3)*DIPE4(3))
      DIPE4D = DIPE4T*AUTODB
      DIPDPT = SQRT(DIPDP(1)*DIPDP(1) + DIPDP(2)*DIPDP(2)
     1              + DIPDP(3)*DIPDP(3))
      DIPDPD = DIPDPT*AUTODB
      WRITE (6,180)
  180 FORMAT (21X,'E4',13X,'DIP',/)
      WRITE (6,190) 'X',DIPE4(1),DIPDP(1)
      WRITE (6,190) 'Y',DIPE4(2),DIPDP(2)
      WRITE (6,190) 'X',DIPE4(3),DIPDP(3)
  190 FORMAT (5X,A1,7X,2F15.6)
      WRITE (6,200) DIPE4T,DIPDPT,
     1               DIPE4D,DIPDPD
  200 FORMAT (//' MAGNITUDE:  ',2F15.6,'  (A.U.)',/,
     1          ' ',12X,2F15.6,'  (DEBYE)')
C
C FIND EIGENVALUES AND EIGENVECTORS OF POLARIZATION MATRIX.
C
      WRITE (6,210)
  210 FORMAT (//,' ',30('*'),' POLARIZABILITY ',20('*'),//)
      WRITE (6,220)
  220 FORMAT (/' E4 POLARIZABILITY TENSOR:')
      I3 = 3
      CALL VECPRT (APOLE4,-I3)
      CALL RSP (APOLE4,I3,I3,EIGS,VECTRS)
      CALL MATOUT (VECTRS,EIGS,I3,I3,I3)
      AVGPE4 = (EIGS(1)+EIGS(2)+EIGS(3))/3.D0
      AVGA3 = AVGPE4*0.14818D00
      AVGESU = AVGPE4*0.296352D-24
      WRITE (6,230)
  230 FORMAT (/' DIP POLARIZABILITY TENSOR:')
      CALL VECPRT (APOLDP,-I3)
      CALL RSP (APOLDP,I3,I3,EIGS,VECTRS)
      CALL MATOUT (VECTRS,EIGS,I3,I3,I3)
      AVGPDP = (EIGS(1)+EIGS(2)+EIGS(3))/3.D0
      AVGA3D = AVGPDP*0.14818D00
      AVGESD = AVGPDP*0.296352D-24
      WRITE (6,240) AVGPE4,AVGPDP,AVGA3,AVGA3D,AVGESU,AVGESD
  240 FORMAT (//,' AVERAGE POLARIZABILITY:',8X,'E4',13X,'DIP',/,
     1           ' ',24X,2F15.6,'  A.U.',/,
     2           ' ',24X,2F15.6,'  ANG.**3',/,
     3           ' ',24X,2(1PD15.6),'  ESU')
C
C  CALCULATE "EXPERIMENTAL" HYPERPOLARIZABILITIES
C
C   8.65710D-33 is a.u. to e.s.u. conversion
      WRITE (6,250)
  250 FORMAT (//,' ',30('*'),' SECOND-ORDER ',25('*'),//)
      BX4 = 0.6D00*(BETAE4(1) + BETAE4(4) + BETAE4(6))
      BY4 = 0.6D00*(BETAE4(2) + BETAE4(5) + BETAE4(8))
      BZ4 = 0.6D00*(BETAE4(3) + BETAE4(7) + BETAE4(9))
      B4MU = (BX4*DIPE4(1) + BY4*DIPE4(2) + BZ4*DIPE4(3))/DIPE4T
      B4ESU = B4MU*8.65710D-03
      BXD = 0.6D00*(BETADP(1) + BETADP(4) + BETADP(6))
      BYD = 0.6D00*(BETADP(2) + BETADP(5) + BETADP(8))
      BZD = 0.6D00*(BETADP(3) + BETADP(7) + BETADP(9))
      BDMU = (BXD*DIPDP(1) + BYD*DIPDP(2) + BZD*DIPDP(3))/DIPDPT
      BDESU = BDMU*8.65710D-03
C
      WRITE (6,260)
  260 FORMAT ('  COMPONENT',12X,'E4',13X,'DIP',/)
      WRITE (6,270) 'XXX',BETAE4(1),BETADP(1)
      WRITE (6,270) 'XYY',BETAE4(4),BETADP(4)
      WRITE (6,270) 'XZZ',BETAE4(6),BETADP(6)
      WRITE (6,270) 'YYY',BETAE4(2),BETADP(2)
      WRITE (6,270) 'YXX',BETAE4(5),BETADP(5)
      WRITE (6,270) 'YZZ',BETAE4(8),BETADP(8)
      WRITE (6,270) 'ZZZ',BETAE4(3),BETADP(3)
      WRITE (6,270) 'ZXX',BETAE4(7),BETADP(7)
      WRITE (6,270) 'ZYY',BETAE4(9),BETADP(9)
  270 FORMAT (' ',5X,A4,5X,2F15.6)
      WRITE (6,280)
  280 FORMAT (/)
      WRITE (6,290) 'BX',BX4,BXD
      WRITE (6,290) 'BY',BY4,BYD
      WRITE (6,290) 'BZ',BZ4,BZD
  290 FORMAT (' ',6X,A2,6X,2F15.6)
      WRITE (6,280)
      WRITE (6,300) B4MU,BDMU,B4ESU,BDESU
  300 FORMAT (' ',4X,'B(AU)',5X,2F15.6,/,
     1        ' ',4X,'B(ESU)',4X,2F15.6,3X,'(X10-30)')
C
      WRITE (6,310)
  310 FORMAT (//' ',30('*'),' THIRD-ORDER ',25('*'),//)
      GAMVAL = (GAMME4(1) + GAMME4(2) + GAMME4(3))
      GAMVAL = GAMVAL + 2.0D00*(GAMME4(4) + GAMME4(5) + GAMME4(6))
      GAMVAL = GAMVAL/5.0D00
C  5.05116D-40 is the a.u. to e.s.u. conversion
      GAMESU = GAMVAL*5.05116D-04
      GAMDIP = (GAMMDP(1) + GAMMDP(2) + GAMMDP(3))
      GAMDIP = GAMDIP + 2.0D00*(GAMMDP(4) + GAMMDP(5) + GAMMDP(6))
      GAMDIP = GAMDIP/5.0D00
      GAMDES = GAMDIP*5.05116D-04
      WRITE (6,320)
  320 FORMAT (' ',17X,'E4',13X,'DIP',/)
      WRITE (6,330) 'XXXX',GAMME4(1),GAMMDP(1)
      WRITE (6,330) 'YYYY',GAMME4(2),GAMMDP(2)
      WRITE (6,330) 'ZZZZ',GAMME4(3),GAMMDP(3)
      WRITE (6,330) 'XXYY',GAMME4(4),GAMMDP(4)
      WRITE (6,330) 'XXZZ',GAMME4(5),GAMMDP(5)
      WRITE (6,330) 'YYZZ',GAMME4(6),GAMMDP(6)
  330 FORMAT (5X,A4,2F15.6)
      WRITE (6,340) GAMVAL,GAMDIP,GAMESU,GAMDES
  340 FORMAT (//' GAMMA = ',1PD15.6,1PD15.6,'  A.U.'/,
     1       ' ',8X,1PD15.6,1PD15.6,'  ESU (X10-36)')
C
      RETURN
      END
      SUBROUTINE FOCK2(F, PTOT, P, W, WJ, WK, NUMAT, NFIRST,
     1NMIDLE, NLAST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION F(*), PTOT(*), WJ(*), WK(*), NFIRST(*), NMIDLE(*),
     1          NLAST(*), P(*), W(*)
C***********************************************************************
C
C THIS IS THE NEW VERSION OF FOCK2 RECEIVED FROM J.J.P.STEWART ON 9JUN89
C 11JUN89 JJPS MOD - "DO 240 I=IA,IC" CHANGED TO "DO 240 I=IA,IB"
C
C FOCK2 FORMS THE TWO-ELECTRON TWO-CENTER REPULSION PART OF THE FOCK
C MATRIX
C ON INPUT  PTOT = TOTAL DENSITY MATRIX.
C           P    = ALPHA OR BETA DENSITY MATRIX.
C           W    = TWO-ELECTRON INTEGRAL MATRIX.
C
C  ON OUTPUT F   = PARTIAL FOCK MATRIX
C***********************************************************************
      COMMON /EULER / TVEC(3,3), ID
      COMMON /KEYWRD/ KEYWRD
      DIMENSION IFACT(MAXORB),
     1I1FACT(MAXORB), JINDEX(256), KINDEX(256), IJPERM(10), LLPERM(10),
     2PK(16), PJA(16), PJB(16), GJA(4,4), GJB(4,4), MMPERM(10),
     3GK(16), IKPERM(16), JPERMA(16), JPERMB(16), P2(MAXORB,MAXORB),
     4F2(MAXORB,MAXORB), PTOT2(MAXORB,MAXORB), JJNDEX(256)
      CHARACTER*80 KEYWRD
      LOGICAL LID
         SAVE                                                           GL0892
      DATA ITYPE /1/
   10 CONTINUE
      GOTO (20,330,50) ITYPE
   20 CONTINUE
C
C   SET UP ARRAY OF LOWER HALF TRIANGLE INDICES (PASCAL'S TRIANGLE)
C
      DO 30 I=1,MAXORB
         IFACT(I)=(I*(I-1))/2
   30 I1FACT(I)=IFACT(I)+I
C
C   SET UP GATHER-SCATTER TYPE ARRAYS FOR USE WITH TWO-ELECTRON
C   INTEGRALS.  JINDEX ARE THE INDICES OF THE J-INTEGRALS FOR ATOM I
C   INTEGRALS.  JJNDEX ARE THE INDICES OF THE J-INTEGRALS FOR ATOM J
C               KINDEX ARE THE INDICES OF THE K-INTEGRALS
C
      M=0
      DO 40 I=1,4
         DO 40 J=1,4
            IJ=MIN(I,J)
            JI=I+J-IJ
            DO 40 K=1,4
               IK=MIN(I,K)
               KI=I+K-IK
               DO 40 L=1,4
                  M=M+1
                  KL=MIN(K,L)
                  LK=K+L-KL
                  JL=MIN(J,L)
                  LJ=J+L-JL
                  KINDEX(M)= IFACT(LJ) +JL + 10*( IFACT(KI) +IK) -10
   40 JINDEX(M)=(IFACT(JI) + IJ)*10 + IFACT(LK) + KL - 10
      L=0
      DO 41 I=1,4
      I1=(I-1)*4
      DO 41 J=1,I
      I1=I1+1
      L=L+1
      IJPERM(L)=I1
      MMPERM(L)=IJPERM(L)-16
      LLPERM(L)=(I1-1)*16
  41  CONTINUE
      L=0
      DO 42 I=1,10
      M=MMPERM(I)
      L=LLPERM(I)
      DO 42 K=1,16
      L=L+1
      M=M+16
  42  JJNDEX(L)=JINDEX(M)
      LID=(ID.EQ.0)
      IONE=1
      IF(ID.NE.0)IONE=0
      IF(INDEX(KEYWRD,'MINDO') .NE. 0) THEN
         ITYPE=2
      ELSE
         ITYPE=3
      ENDIF
      GOTO 10
   50 KK=0
      L=0
      DO 51 I=1,NUMAT
      IA=NFIRST(I)
      IB=NLAST(I)
      M=0
      DO 141 J=IA,IB
      DO 141 K=IA,IB
      M=M+1
      JK=MIN(J,K)
      KJ=K+J-JK
      JK=JK+(KJ*(KJ-1))/2
      PTOT2(I,M)=PTOT(JK)
 141  CONTINUE
  51  CONTINUE
      DO 320 II=1,NUMAT
         IA=NFIRST(II)
         IB=NLAST(II)
         IMINUS=II-IONE
         DO 310 JJ=1,IMINUS
            JA=NFIRST(JJ)
            JB=NLAST(JJ)
            JC=NMIDLE(JJ)
            IF(LID) THEN
               DREP=W(KK+1)
               IF(IB-IA.GE.3.AND.JB-JA.GE.3)THEN
C
C                         HEAVY-ATOM  - HEAVY-ATOM
C
C   EXTRACT COULOMB TERMS
C
                  DO 311 I=1,16
                  PJA(I)=PTOT2(II,I)
  311             PJB(I)=PTOT2(JJ,I)
C
C  COULOMB TERMS
C
               CALL JAB(IA,JA,LLPERM,JINDEX, JJNDEX, PJA,PJB,W(KK+1),F)
C
C  EXCHANGE TERMS
C
C
C  EXTRACT INTERSECTION OF ATOMS II AND JJ IN THE SPIN DENSITY MATRIX
C
                  L=0
                  DO 110 I=IA,IB
                     I1=IFACT(I)+JA
                     DO 110 J=I1,I1+3
                        L=L+1
  110             PK(L)=P(J)
                  CALL KAB(IA,JA, PK, W(KK+1), KINDEX, F)
                  KK=KK+100
               ELSEIF(IB-IA.GE.3)THEN
C
C                         LIGHT-ATOM  - HEAVY-ATOM
C
C
C   COULOMB TERMS
C
                  SUMDIA=0.D0
                  SUMOFF=0.D0
                  LL=I1FACT(JA)
                  K=0
                  DO 211 I=0,3
                     J1=IFACT(IA+I)+IA-1
                     DO 201 J=0,I-1
                        K=K+1
                        J1=J1+1
                        F(J1)=F(J1)+PTOT(LL)*W(KK+K)
  201                SUMOFF=SUMOFF+PTOT(J1)*W(KK+K)
                     J1=J1+1
                     K=K+1
                     F(J1)=F(J1)+PTOT(LL)*W(KK+K)
  211             SUMDIA=SUMDIA+PTOT(J1)*W(KK+K)
                  F(LL)=F(LL)+SUMOFF*2.D0+SUMDIA
C
C  EXCHANGE TERMS
C
C
C  EXTRACT INTERSECTION OF ATOMS II AND JJ IN THE SPIN DENSITY MATRIX
C
                  K=0
                  DO 231 I=IA,IB
                  I1=IFACT(I)+JA
                     SUM=0.D0
                     DO 221 J=IA,IB
                     K=K+1
                        J1=IFACT(J)+JA
  221                SUM=SUM+P(J1)*W(KK+JINDEX(K))
  231             F(I1)=F(I1)-SUM
                  KK=KK+10
               ELSEIF(JB-JA.GE.3)THEN
C
C                         HEAVY-ATOM - LIGHT-ATOM
C
C
C   COULOMB TERMS
C
                  SUMDIA=0.D0
                  SUMOFF=0.D0
                  LL=I1FACT(IA)
                  K=0
                  DO 210 I=0,3
                     J1=IFACT(JA+I)+JA-1
                     DO 200 J=0,I-1
                        K=K+1
                        J1=J1+1
                        F(J1)=F(J1)+PTOT(LL)*W(KK+K)
  200                SUMOFF=SUMOFF+PTOT(J1)*W(KK+K)
                     J1=J1+1
                     K=K+1
                     F(J1)=F(J1)+PTOT(LL)*W(KK+K)
  210             SUMDIA=SUMDIA+PTOT(J1)*W(KK+K)
                  F(LL)=F(LL)+SUMOFF*2.D0+SUMDIA
C
C  EXCHANGE TERMS
C
C
C  EXTRACT INTERSECTION OF ATOMS II AND JJ IN THE SPIN DENSITY MATRIX
C
                  K=IFACT(IA)+JA
                  J=0
                  DO 230 I=K,K+3
                     SUM=0.D0
                     DO 220 L=K,K+3
                        J=J+1
  220                SUM=SUM+P(L)*W(KK+JINDEX(J))
  230             F(I)=F(I)-SUM
                  KK=KK+10
               ELSE
C
C                         LIGHT-ATOM - LIGHT-ATOM
C
               I1=I1FACT(IA)
               J1=I1FACT(JA)
               IJ=I1+JA-IA
               F(I1)=F(I1)+PTOT(J1)*W(KK+1)
               F(J1)=F(J1)+PTOT(I1)*W(KK+1)
               F(IJ)=F(IJ)-P   (IJ)*W(KK+1)
               KK=KK+1
               ENDIF
            ELSE
               DREP=WJ(KK+1)
               DO 240 I=IA,IB
                  KA=IFACT(I)
                  DO 240 J=IA,I
                     KB=IFACT(J)
                     IJ=KA+J
                     AA=2.0D00
                     IF (I.EQ.J) AA=1.0D00
                     DO 240 K=JA,JC
                        KC=IFACT(K)
                        IF(I.GE.K) THEN
                           IK=KA+K
                        ELSE
                           IK=0
                        ENDIF
                        IF(J.GE.K) THEN
                           JK=KB+K
                        ELSE
                           JK=0
                        ENDIF
                        DO 240 L=JA,K
                           IF(I.GE.L) THEN
                              IL=KA+L
                           ELSE
                              IL=0
                           ENDIF
                           IF(J.GE.L) THEN
                              JL=KB+L
                           ELSE
                              JL=0
                           ENDIF
                           KL=KC+L
                           BB=2.0D00
                           IF (K.EQ.L) BB=1.0D00
                           KK=KK+1
                           AJ=WJ(KK)
                           AK=WK(KK)
C
C     A  IS THE REPULSION INTEGRAL (I,J/K,L) WHERE ORBITALS I AND J ARE
C     ON ATOM II, AND ORBITALS K AND L ARE ON ATOM JJ.
C     AA AND BB ARE CORRECTION FACTORS SINCE
C     (I,J/K,L)=(J,I/K,L)=(I,J/L,K)=(J,I/L,K)
C     IJ IS THE LOCATION OF THE MATRIX ELEMENTS BETWEEN ATOMIC ORBITALS
C     I AND J.  SIMILARLY FOR IK ETC.
C
C THIS FORMS THE TWO-ELECTRON TWO-CENTER REPULSION PART OF THE FOCK
C MATRIX.  THE CODE HERE IS HARD TO FOLLOW, AND IMPOSSIBLE TO MODIFY!,
C BUT IT WORKS,
                           IF(KL.LE.IJ)THEN
                              IF(I.EQ.K.AND.AA+BB.LT.2.1D0)THEN
                                 BB=BB*0.5D0
                                 AA=AA*0.5D0
                                 F(IJ)=F(IJ)+BB*AJ*PTOT(KL)
                                 F(KL)=F(KL)+AA*AJ*PTOT(IJ)
                              ELSE
                                 F(IJ)=F(IJ)+BB*AJ*PTOT(KL)
                                 F(KL)=F(KL)+AA*AJ*PTOT(IJ)
                                 A=AK*AA*BB*0.25D0
                                 F(IK)=F(IK)-A*P(JL)
                                 F(IL)=F(IL)-A*P(JK)
                                 F(JK)=F(JK)-A*P(IL)
                                 F(JL)=F(JL)-A*P(IK)
                              ENDIF
                           ENDIF
  240          CONTINUE
            ENDIF
  310    CONTINUE
  320 CONTINUE
C
      RETURN
  330 KR=0
      DO 360 II=1,NUMAT
         IA=NFIRST(II)
         IB=NLAST(II)
         IM1=II-IONE
         DO 350 JJ=1,IM1
            KR=KR+1
            IF(LID)THEN
               ELREP=W(KR)
               ELEXC=ELREP
            ELSE
               ELREP=WJ(KR)
               ELEXC=WK(KR)
            ENDIF
            JA=NFIRST(JJ)
            JB=NLAST(JJ)
            DO 340 I=IA,IB
               KA=IFACT(I)
               KK=KA+I
               DO 340 K=JA,JB
                  LL=I1FACT(K)
                  IK=KA+K
                  F(KK)=F(KK)+PTOT(LL)*ELREP
                  F(LL)=F(LL)+PTOT(KK)*ELREP
  340       F(IK)=F(IK)-P(IK)*ELEXC
  350    CONTINUE
  360 CONTINUE
      RETURN
      END
      SUBROUTINE FOCK2D(F,PTOT, P, W, WJ, WK, NUMAT, NFIRST,
     1NMIDLE, NLAST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(*), PTOT(*), WJ(*), WK(*), NFIRST(*), NMIDLE(*),
     1          NLAST(*), P(*), W(*)
C***********************************************************************
C
C FOCK2 FORMS THE TWO-ELECTRON TWO-CENTER REPULSION PART OF THE FOCK
C MATRIX
C ON INPUT  PTOT = TOTAL DENSITY MATRIX.
C           P    = ALPHA OR BETA DENSITY MATRIX.
C           W    = TWO-ELECTRON INTEGRAL MATRIX.
C
C  ON OUTPUT F   = PARTIAL FOCK MATRIX
C***********************************************************************
      COMMON /EULER / TVEC(3,3), ID
      COMMON /KEYWRD/ KEYWRD
      DIMENSION SPPOP(2), DPOP(2), IFACT(20), I1FACT(20)
      LOGICAL LID
      CHARACTER*80 KEYWRD
         SAVE                                                           GL0892
      DATA ITYPE /1/
   10 CONTINUE
      GOTO (20,170,40) ITYPE
   20 CONTINUE
C
C   SET UP ARRAY OF (I*(I-1))/2
C
      DO 30 I=1,20
         IFACT(I)=(I*(I-1))/2
   30 I1FACT(I)=IFACT(I)+I
      LID=(ID.EQ.0)
      IONE=1
      IF(INDEX(KEYWRD,'MINDO') .NE. 0) THEN
         ITYPE=2
      ELSE
         ITYPE=3
      ENDIF
      GOTO 10
   40 KK=0
      NORBS=NLAST(NUMAT)
      LINEA1=(NORBS*(NORBS+1))/2 + 1
      P(LINEA1)=0.D0
      DO 160 II=1,NUMAT
         IA=NFIRST(II)
         IB=NLAST(II)
         IC=NMIDLE(II)
         SUM=0.D0
         DO 50 I=IA,IC
   50    SUM=SUM+PTOT(I1FACT(I))
         SPPOP(II)=SUM
         SUM=0.D0
         DO 60 I=IC+1,IB
   60    SUM=SUM+PTOT(I1FACT(I))
         DPOP(II)=SUM
         IMINUS=II-IONE
         DO 150 JJ=1,IMINUS
            JA=NFIRST(JJ)
            JB=NLAST(JJ)
            JC=NMIDLE(JJ)
            DREP=WJ(KK+1)
            IF(LID) THEN
               DREP=W(KK+1)
               DO 70 I=IA,IC
                  KA=IFACT(I)
                  DO 70 J=IA,I
                     KB=IFACT(J)
                     IJ=KA+J
                     AA=2.0D00
                     IF (I.EQ.J) AA=1.0D00
                     DO 70 K=JA,JC
                        KC=IFACT(K)
                        IK=KA+K
                        JK=KB+K
                        DO 70 L=JA,K
                           IL=KA+L
                           JL=KB+L
                           KL=KC+L
                           BB=2.0D00
                           IF (K.EQ.L) BB=1.0D00
                           KK=KK+1
                           A=W(KK)
C
C     A  IS THE REPULSION INTEGRAL (I,J/K,L) WHERE ORBITALS I AND J ARE
C     ON ATOM II, AND ORBITALS K AND L ARE ON ATOM JJ.
C     AA AND BB ARE CORRECTION FACTORS SINCE
C     (I,J/K,L)=(J,I/K,L)=(I,J/L,K)=(J,I/L,K)
C     IJ IS THE LOCATION OF THE MATRIX ELEMENTS BETWEEN ATOMIC ORBITALS
C     I AND J.  SIMILARLY FOR IK ETC.
C
C THIS FORMS THE TWO-ELECTRON TWO-CENTER REPULSION PART OF THE FOCK
C MATRIX.  THE CODE HERE IS HARD TO FOLLOW, AND IMPOSSIBLE TO MODIFY!,
C BUT IT WORKS,
                           F(IJ)=F(IJ)+BB*A*PTOT(KL)
                           F(KL)=F(KL)+AA*A*PTOT(IJ)
                           A=A*AA*BB*0.25D0
                           F(IK)=F(IK)-A*P(JL)
                           F(IL)=F(IL)-A*P(JK)
                           F(JK)=F(JK)-A*P(IL)
                           F(JL)=F(JL)-A*P(IK)
   70          CONTINUE
            ELSE
               DO 80 I=IA,IC
                  KA=IFACT(I)
                  DO 80 J=IA,I
                     KB=IFACT(J)
                     IJ=KA+J
                     AA=2.0D00
                     IF (I.EQ.J) AA=1.0D00
                     DO 80 K=JA,JC
                        KC=IFACT(K)
                        IF(I.GE.K) THEN
                           IK=KA+K
                        ELSE
                           IK=LINEA1
                        ENDIF
                        IF(J.GE.K) THEN
                           JK=KB+K
                        ELSE
                           JK=LINEA1
                        ENDIF
                        DO 80 L=JA,K
                           IF(I.GE.L) THEN
                              IL=KA+L
                           ELSE
                              IL=LINEA1
                           ENDIF
                           IF(J.GE.L) THEN
                              JL=KB+L
                           ELSE
                              JL=LINEA1
                           ENDIF
                           KL=KC+L
                           BB=2.0D00
                           IF (K.EQ.L) BB=1.0D00
                           KK=KK+1
                           AJ=WJ(KK)
                           AK=WK(KK)
C
C     A  IS THE REPULSION INTEGRAL (I,J/K,L) WHERE ORBITALS I AND J ARE
C     ON ATOM II, AND ORBITALS K AND L ARE ON ATOM JJ.
C     AA AND BB ARE CORRECTION FACTORS SINCE
C     (I,J/K,L)=(J,I/K,L)=(I,J/L,K)=(J,I/L,K)
C     IJ IS THE LOCATION OF THE MATRIX ELEMENTS BETWEEN ATOMIC ORBITALS
C     I AND J.  SIMILARLY FOR IK ETC.
C
C THIS FORMS THE TWO-ELECTRON TWO-CENTER REPULSION PART OF THE FOCK
C MATRIX.  THE CODE HERE IS HARD TO FOLLOW, AND IMPOSSIBLE TO MODIFY!,
C BUT IT WORKS,
                           IF(KL.LE.IJ)THEN
                              IF(I.EQ.K.AND.AA+BB.LT.2.1D0)THEN
                                 BB=BB*0.5D0
                                 AA=AA*0.5D0
                                 F(IJ)=F(IJ)+BB*AJ*PTOT(KL)
                                 F(KL)=F(KL)+AA*AJ*PTOT(IJ)
                              ELSE
                                 F(IJ)=F(IJ)+BB*AJ*PTOT(KL)
                                 F(KL)=F(KL)+AA*AJ*PTOT(IJ)
                                 A=AK*AA*BB*0.25D0
                                 F(IK)=F(IK)-A*P(JL)
                                 F(IL)=F(IL)-A*P(JK)
                                 F(JK)=F(JK)-A*P(IL)
                                 F(JL)=F(JL)-A*P(IK)
                              ENDIF
                           ENDIF
   80          CONTINUE
            ENDIF
C
C   D-ORBITAL CORRECTION
C
            DO 90 I=IC+1,IB
               KA=IFACT(I)
               DO 90 J=JA,JB
                  IJ=KA+J
C
C   ATOM J (S, P, AND D (IF PRESENT)) EXCHANGE WITH ATOM I (D ONLY)
C
   90       F(IJ)=F(IJ)-0.5D0*DREP*P(IJ)
            DO 100 I=IA,IC
               KA=IFACT(I)
               DO 100 J=JC+1,JB
                  IJ=KA+J
C
C    ATOM J (D(IF PRESENT)) EXCHANGE WITH ATOM I (S AND P ONLY)
C
  100       F(IJ)=F(IJ)-0.5D0*DREP*P(IJ)
C
C                      THE COULOMB REPULSION TERMS.
C
C     FIRST, ATOM J (S, P AND D SHELLS) BEING REPELLED BY ATOM I(DSHELL)
C
cdir$ ivdep
            DO 110 J=JA,JB
               J2=I1FACT(J)
  110       F(J2)=F(J2)+DREP*DPOP(II)
C
C     ATOM J (D SHELL) BEING REPELLED BY ATOM I (S AND P SHELLS)
C
cdir$ ivdep
            DO 120 J=JC+1,JB
               J2=I1FACT(J)
  120       F(J2)=F(J2)+DREP*SPPOP(II)
C
C     ATOM I (S, P AND D SHELLS) BEING REPELLED BY ATOM J (D SHELL)
C
cdir$ ivdep
            DO 130 I=IA,IB
               I2=I1FACT(I)
  130       F(I2)=F(I2)+DREP*DPOP(JJ)
C
C    ATOM I (D SHELL) BEING REPELLED BY ATOM J (S AND P SHELLS)
C
cdir$ ivdep
            DO 140 I=IC+1,IB
               I2=I1FACT(I)
  140       F(I2)=F(I2)+DREP*SPPOP(JJ)
  150    CONTINUE
  160 CONTINUE
C
      RETURN
  170 KR=0
      DO 200 II=1,NUMAT
         IA=NFIRST(II)
         IB=NLAST(II)
         IM1=II-IONE
         DO 190 JJ=1,IM1
            KR=KR+1
            IF(LID)THEN
               ELREP=W(KR)
               ELEXC=ELREP
            ELSE
               ELREP=WJ(KR)
               ELEXC=WK(KR)
            ENDIF
            JA=NFIRST(JJ)
            JB=NLAST(JJ)
            DO 180 I=IA,IB
               KA=IFACT(I)
               KK=KA+I
               DO 180 K=JA,JB
                  LL=I1FACT(K)
                  IK=KA+K
                  F(KK)=F(KK)+PTOT(LL)*ELREP
                  F(LL)=F(LL)+PTOT(KK)*ELREP
  180       F(IK)=F(IK)-P(IK)*ELEXC
  190    CONTINUE
  200 CONTINUE
      RETURN
      END
      SUBROUTINE FORSAV(TIME,DELDIP,IPT,N3,FMATRX, COORD,NVAR,REFH,
     1                  EVECS,JSTART,FCONST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION FMATRX(*), DELDIP(3,*), COORD(*), EVECS(*), FCONST(*)
************************************************************************
*
*  FORSAV SAVES AND RESTORES DATA USED IN THE FORCE CALCULATION.
*
* ON INPUT TIME = TOTAL TIME ELAPSED SINCE THE START OF THE CALCULATION.
*          IPT  = LINE OF FORCE MATRIX REACHED, IF IN WRITE MODE,
*               = 0 IF IN READ MODE.
*        FMATRX = FORCE MATRIX
************************************************************************
      COMMON /DENSTY/ P(MPACK), PA(MPACK), PB(MPACK)
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
         SAVE                                                           GL0892
      OPEN(UNIT=9,FILE='FOR009',STATUS='UNKNOWN',FORM='UNFORMATTED')
      REWIND 9
      OPEN(UNIT=10,FILE='FOR010',STATUS='UNKNOWN',FORM='UNFORMATTED')
      REWIND 10
      IR=9
      IW=9
      IF( IPT .EQ. 0 ) THEN
C
C   READ IN FORCE DATA
C
         READ(IR,END=20,ERR=20)TIME,IPT,REFH
         LINEAR=(NVAR*(NVAR+1))/2
         READ(IR)(COORD(I),I=1,NVAR)
         READ(IR,END=10,ERR=10)(FMATRX(I),I=1,LINEAR)
         READ(IR)((DELDIP(J,I),J=1,3),I=1,IPT)
         N33=NVAR*NVAR
         READ(IR)(EVECS(I),I=1,N33)
         READ(IR)JSTART,(FCONST(I),I=1,NVAR)
         RETURN
      ELSE
C
C    WRITE FORCE DATA
C
         REWIND IW
         IF(TIME.GT.1.D6)TIME=TIME-1.D6
         WRITE(IW)TIME,IPT,REFH
         LINEAR=(NVAR*(NVAR+1))/2
         WRITE(IW)(COORD(I),I=1,NVAR)
         WRITE(IW)(FMATRX(I),I=1,LINEAR)
         WRITE(IW)((DELDIP(J,I),J=1,3),I=1,IPT)
         N33=NVAR*NVAR
         WRITE(IR)(EVECS(I),I=1,N33)
         WRITE(IR)JSTART,(FCONST(I),I=1,NVAR)
         LINEAR=(NORBS*(NORBS+1))/2
         WRITE(10)(PA(I),I=1,LINEAR)
         IF(NALPHA.NE.0)WRITE(10)(PB(I),I=1,LINEAR)
         CLOSE(9)
         CLOSE(10)
      ENDIF
      RETURN
   10 WRITE(6,'(''INSUFFICIENT DATA ON DISK FILES FOR A FORCE '',
     1''CALCULATION'',/10X,''RESTART. PERHAPS THIS STARTED OF AS A '',
     2''FORCE CALCULATION '')')
      WRITE(6,'(10X,''BUT THE GEOMETRY HAD TO BE OPTIMIZED FIRST, '',
     1''IN WHICH CASE '',/10X,''REMOVE THE KEY-WORD "FORCE".'')')
      STOP
   20 WRITE(6,'(//10X,''NO RESTART FILE EXISTS!'')')
      STOP
      END
      SUBROUTINE FRAME(FMAT,NUMAT,MODE,SHIFT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION FMAT(*), SHIFT(6)
      COMMON /COORD /COORD(3,NUMATM)
      COMMON /ATMASS/ ATMASS(NUMATM)
      DIMENSION VIB(6,MAXPAR), ROT(3,3), COORD1(3,NUMATM)
***********************************************************************
*
*   FRAME APPLIES AN RIGID ORIENTATION TO THE MOLECULE IN A FORCE
*         CALCULATION. THE TRANSLATIONS ARE GIVEN A 'FORCE CONSTANT'
*         OF T(X)=500 MILLIDYNES/ANGSTROM
*            T(Y)=600 MILLIDYNES/ANGSTROM
*            T(Z)=700 MILLIDYNES/ANGSTROM
*         AND THE ROTATIONS ARE GIVEN A 'FORCE CONSTANT' OF
*            R(X)=800 MILLIDYNES/ANGSTROM
*            R(Y)=900 MILLIDYNES/ANGSTROM
*            R(Z)=1000 MILLIDYNES/ANGSTROM,
*    THE ROTATIONS ARE MADE ABOUT AXES DETERMINED BY THE MOMENTS
*    OF INERTIA, WHICH IN TURN DEPEND ON THE ISOTOPIC MASSES. FOR
*    THE NORMAL FREQUENCY CALCULATION THESE ARE THE REAL MASSES,
*    FOR THE FORCE CALCULATION THEY ARE ALL UNITY.
***********************************************************************
      COMMON /EULER / TVEC(3,3), ID
         SAVE                                                           GL0892
      CALL AXIS(COORD,NUMAT,A,B,C,SUMW, MODE,ROT )
      DO 20 I=1,NUMAT
         DO 20 J=1,3
            SUM=0.D0
            DO 10 K=1,3
   10       SUM=SUM+COORD(K,I)*ROT(K,J)
   20 COORD1(J,I)=SUM
      N3=NUMAT*3
      J=0
      WTMASS=1.D0
      DO 30 I=1,NUMAT
         IF(MODE.EQ.1)  WTMASS=SQRT(ATMASS(I))
         J=J+1
         VIB(1,J)=WTMASS
         VIB(2,J)=0.D0
         VIB(3,J)=0.D0
         VIB(4,J)=0.D0
         VIB(5,J)=COORD1(3,I)*WTMASS
         VIB(6,J)=COORD1(2,I)*WTMASS
         J=J+1
         VIB(1,J)=0.D0
         VIB(2,J)=WTMASS
         VIB(3,J)=0.D0
         VIB(4,J)=COORD1(3,I)*WTMASS
         VIB(5,J)=0.D0
         VIB(6,J)=-COORD1(1,I)*WTMASS
         J=J+1
         VIB(1,J)=0.D0
         VIB(2,J)=0.D0
         VIB(3,J)=WTMASS
         VIB(4,J)=-COORD1(2,I)*WTMASS
         VIB(5,J)=-COORD1(1,I)*WTMASS
         VIB(6,J)=0.D0
   30 CONTINUE
      J=1
      DO 50 I=1,NUMAT
         DO 40 K=4,6
            X=VIB(K,J)
            Y=VIB(K,J+1)
            Z=VIB(K,J+2)
            VIB(K,J  )=X*ROT(1,1)+Y*ROT(1,2)+Z*ROT(1,3)
            VIB(K,J+1)=X*ROT(2,1)+Y*ROT(2,2)+Z*ROT(2,3)
            VIB(K,J+2)=X*ROT(3,1)+Y*ROT(3,2)+Z*ROT(3,3)
   40    CONTINUE
         J=J+3
   50 CONTINUE
      SUM1=0.D0
      SUM2=0.D0
      SUM3=0.D0
      SUM4=0.D0
      SUM5=0.D0
      SUM6=0.D0
      DO 60 I=1,N3
         SUM1=SUM1+VIB(1,I)**2
         SUM2=SUM2+VIB(2,I)**2
         SUM3=SUM3+VIB(3,I)**2
         SUM4=SUM4+VIB(4,I)**2
         SUM5=SUM5+VIB(5,I)**2
   60 SUM6=SUM6+VIB(6,I)**2
      IF(SUM1.GT.1.D-5)SUM1=SQRT(1.D0/SUM1)
      IF(SUM2.GT.1.D-5)SUM2=SQRT(1.D0/SUM2)
      IF(SUM3.GT.1.D-5)SUM3=SQRT(1.D0/SUM3)
      IF(SUM4.GT.1.D-5)SUM4=SQRT(1.D0/SUM4)
      IF(SUM5.GT.1.D-5)SUM5=SQRT(1.D0/SUM5)
      IF(SUM6.GT.1.D-5)SUM6=SQRT(1.D0/SUM6)
      IF(ID.NE.0)THEN
         SUM4=0.D0
         SUM5=0.D0
         SUM6=0.D0
      ENDIF
      DO 70 I=1,N3
         VIB(1,I)=VIB(1,I)*SUM1
         VIB(2,I)=VIB(2,I)*SUM2
         VIB(3,I)=VIB(3,I)*SUM3
         VIB(4,I)=VIB(4,I)*SUM4
         VIB(5,I)=VIB(5,I)*SUM5
   70 VIB(6,I)=VIB(6,I)*SUM6
      DO 80 I=1,6
   80 SHIFT(I)=400.D0+I*100.D0
      L=0
      DO 100 I=1,N3
         DO 100 J=1,I
            L=L+1
            SUM1=0.D0
            DO 90 K=1,6
   90       SUM1=SUM1+VIB(K,I)*SHIFT(K)*VIB(K,J)
  100 FMAT(L)=FMAT(L)+SUM1
      END
      SUBROUTINE FREQCY(FMATRX,FREQ,CNORML,REDMAS,TRAVEL,EORC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION FMATRX(*), REDMAS(*), FREQ(*), CNORML(*), TRAVEL(*)
      LOGICAL EORC
*********************************************************************
*
*  FRCE CALCULATES THE FORCE CONSTANTS AND VIBRATIONAL FREQUENCIES
*       FOR A MOLECULE.  IT USES THE ISOTOPIC MASSES TO WEIGHT THE
*       FORCE MATRIX
*
* ON INPUT   FMATRX   =  FORCE MATRIX, OF SIZE NUMAT*3*(NUMAT*3+1)/2.
*
*********************************************************************
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /ATMASS/ ATMASS(NUMATM)
      COMMON /SCRACH/ OLDF(MAXHES), DUMMY(MAXPAR**2-MAXHES)
      DIMENSION WTMASS(MAXPAR), SHIFT(6)
         SAVE                                                           GL0892
      DATA FACT/6.023D23/
C
C    CONVERSION FACTOR FOR SPEED OF LIGHT AND 2 PI.
C
      C2PI=1.D0/(2.998D10*3.141592653598D0*2.D0)
C NOW TO CALCULATE THE VIBRATIONAL FREQUENCIES
C
C   FIND CONVERSION CONSTANTS FOR MASS WEIGHTED SYSTEM
      N3=3*NUMAT
      L=0
      DO 10 I=1,NUMAT
         WEIGHT=1.4142136D0/SQRT(ATMASS(I))
         DO 10 J=1,3
            L=L+1
   10 WTMASS(L)=WEIGHT
C    CONVERT TO MASS WEIGHTED FMATRX
      LINEAR=0
      DO 20 I=1,N3
         DO 20 J=1,I
            LINEAR=LINEAR+1
            OLDF(LINEAR)=  FMATRX(LINEAR)*1.D5
   20 FMATRX(LINEAR)=FMATRX(LINEAR)*WTMASS(I)*WTMASS(J)
C
C    1.D5 IS TO CONVERT FROM MILLIDYNES/ANGSTROM TO DYNES/CM.
C
C    DIAGONALIZE
      CALL FRAME(FMATRX,NUMAT,1, SHIFT)
      CALL RSP(FMATRX,N3,N3,FREQ,CNORML)
      DO 30 I=1,N3
         J=(FREQ(I)+50.D0)*0.01D0
   30 FREQ(I)=FREQ(I)-J*100
      DO 40 I=1,N3
   40 FREQ(I)=FREQ(I)*1.D5
C
C    CALCULATE REDUCED MASSES, STORE IN REDMAS
C
      DO 80 I=1,N3
         II=(I-1)*N3
         SUM=0.D0
         DO 70 J=1,N3
            JII=J+II
            JJ=(J*(J-1))/2
            DO 50 K=1,J
   50       SUM=SUM+CNORML(JII)*OLDF(JJ+K)*CNORML(K+II)
            DO 60 K=J+1,N3
   60       SUM=SUM+CNORML(JII)*OLDF((K*(K-1))/2+J)*CNORML(K+II)
   70    CONTINUE
         SUM1=SUM*2.D0
         IF(ABS(FREQ(I)).GT.ABS(SUM)*1.D-20) THEN
            SUM=1.D0*SUM/FREQ(I)
         ELSE
            SUM=0.D0
         ENDIF
         FREQ(I)=SIGN(SQRT(FACT*ABS(FREQ(I)))*C2PI,FREQ(I))
         IF(ABS(FREQ(I)).LT.ABS(SUM1)*1.D+20) THEN
            SUM1=SQRT(ABS(FREQ(I)/(SUM1*1.D-5)))
         ELSE
            SUM1=0.D0
         ENDIF
         IF(SUM.LT.0.D0.OR.SUM.GT.100)SUM=0.D0
C
C 0.0063024=SQRT(2*A*B*C/N) WHERE
C         A=1.196D8 = CONVERSION OF CM**(-1) TO (ERGS = DYNE.ANGSTROMS)
C         B=1000.0  = MILLIDYNES TO DYNES
C         C=1.D8    = CENTIMETERS TO ANGSTROMS
C         N=6.02205D23 = AVOGADRO'S NUMBER
         TRAVEL(I)=SUM1*0.0063024D0
         IF(TRAVEL(I).GT.1.D0)TRAVEL(I)=0.D0
C#      WRITE(6,*)TRAVEL(I)
   80 REDMAS(I)=SUM
      IF(EORC) THEN
C
C    SWITCH EIGENVALUES TO FREQUENCIES
C
C    CONVERT NORMAL VECTORS TO CARTESIAN COORDINATES
C    AND NORMALIZE SO THAT THE TOTAL MOVEMENT IS 1.0 ANGSTROM.
C
         IJ=0
         DO 120 I=1,N3
            SUM=0.D0
            J=0
            DO 100 JJ=1,NUMAT
               SUM1=0.D0
               DO 90 J1=1,3
                  J=J+1
                  IJ=IJ+1
                  CNORML(IJ)=CNORML(IJ)*WTMASS(J)
   90          SUM1=SUM1+CNORML(IJ)**2
  100       SUM=SUM+SQRT(SUM1)
            SUM=1.D0/SUM
            IJ=IJ-N3
            DO 110 J=1,N3
               IJ=IJ+1
  110       CNORML(IJ)=CNORML(IJ)*SUM
  120    CONTINUE
C
C          RETURN HESSIAN IN MILLIDYNES/ANGSTROM IN FMATRX
C
         DO 130 I=1,LINEAR
  130    FMATRX(I)=OLDF(I)*1.D-5
      ELSE
C
C  RETURN HESSIAN AS MASS-WEIGHTED FMATRIX
         LINEAR=0
C
         DO 140 I=1,N3
            DO 140 J=1,I
               LINEAR=LINEAR+1
  140    FMATRX(LINEAR)=OLDF(LINEAR)*1.D-5*WTMASS(I)*WTMASS(J)
      ENDIF
      RETURN
      END
      SUBROUTINE GEOUT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
**********************************************************************
*
*   GEOUT PRINTS THE CURRENT GEOMETRY.  IT CAN BE CALLED ANY TIME,
*         FROM ANY POINT IN THE PROGRAM AND DOES NOT AFFECT ANYTHING.
*
**********************************************************************
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /GEOVAR/ XPARAM(MAXPAR), NVAR, LOC(2,MAXPAR)               IR0394
C    changed common path for portability  (IR)
      COMMON /PATHI / LATOM,LPARAM
      COMMON /PATHR / REACT(200)
c
      COMMON /ELEMTS/ ELEMNT(107)
      DIMENSION COORD(3,NUMATM)
      CHARACTER Q(3)*1, ELEMNT*2
      LOGICAL CART
         SAVE                                                           GL0892
C
C *** OUTPUT THE PARAMETER DATA.
C
      CART=.FALSE.
      IF(NA(1).NE.0) THEN
         CART=.TRUE.
         CALL XYZINT(GEO,NATOMS,NA,NB,NC,1.D0,COORD)
         NA(1)=99
      ELSE
         DO 10 I=1,NATOMS
            DO 10 J=1,3
   10    COORD(J,I)=GEO(J,I)
      ENDIF
      DEGREE=57.29577951D00
      WRITE (6,20)
   20 FORMAT (/4X,'ATOM',3X,'CHEMICAL',2X,'BOND LENGTH',4X,'BOND ANGLE'
     1,4X ,' TWIST ANGLE',/3X,'NUMBER',2X,'SYMBOL', 4X,'(ANGSTROMS)',5
     2X,'(DEGREES)',5X,' (DEGREES)',/4X,'(I)',18X,'NA:I',10X,'NB:NA:I',5
     3X,' NC:NB:NA:I',5X,'NA',3X,'NB',3X,'NC',/)
      N=1
      IA=LOC(1,1)
      DO 50 I=1,NATOMS
         DO 30 J=1,3
            Q(J)=' '
            IF (IA.NE.I) GO TO 30
            IF (J.NE.LOC(2,N)) GO TO 30
            Q(J)='*'
            N=N+1
            IA=LOC(1,N)
   30    CONTINUE
         W = COORD(2,I) * DEGREE
         X = COORD(3,I) * DEGREE
         IF (LATOM.NE.I) GO TO 40
         J=LPARAM
         Q(J)='+ '
   40    CONTINUE
         IF(LABELS(I).NE.0)THEN
            IF(I.GT.3)THEN
               WRITE (6,'(3X,I4 ,5X,A2,F16.5,1X,A1,F15.5,1X,A1,F12.5,1X,
     1A1,I5,2I5)') I,ELEMNT(LABELS(I)),COORD(1,I),Q(1),W,Q(2),X,Q(3),NA(
     2I),NB(I),NC(I)
            ELSEIF(I.EQ.3)THEN
               WRITE (6,'(''      3'',5X,A2,F16.5,1X,A1,F15.5,1X,A1,14X,
     12I5)') ELEMNT(LABELS(3)),COORD(1,3),Q(1),W,Q(2),NA(3),NB(3)
            ELSEIF(I.EQ.2)THEN
               WRITE (6,'(''      2'',5X,A2,F16.5,1X,A1,31X,I5)')
     1 ELEMNT(LABELS(2)),COORD(1,2),Q(1),NA(2)
            ELSE
               WRITE (6,'(''      1'',5X,A2)') ELEMNT(LABELS(1))
            ENDIF
         ENDIF
   50 CONTINUE
      RETURN
      END
      SUBROUTINE GOVER(NI,NJ,XI,XJ,R,SG)
************************************************************************
*                                                                      *
*   GOVER CALCULATES THE OVERLAP INTEGRALS USING A GAUSSIAN EXPANSION  *
*         STO-6G BY R.F. STEWART, J. CHEM. PHYS., 52 431-438, 1970     *
*                                                                      *
*         ON INPUT   NI   =  ATOMIC NUMBER OF FIRST ATOM               *
*                    NJ   =  ATOMIC NUMBER OF SECOND ATOM              *
*                    R    =  INTERATOMIC DISTANCE IN ANGSTROMS         *
*         ON EXIT    S    =  9X9 ARRAY OF OVERLAPS, IN ORDER S,PX,PY,  *
*                            PZ                                        *
*                                                                      *
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      COMMON /NATYPE/ NZTYPE(107), MTYPE(30),LTYPE
      COMMON /TEMP/  C(60,6), Z(60,6)
      COMMON /NATORB/ NATORB(107)
      DIMENSION S(6,6), XI(3),XJ(3), SG(9,9)
         SAVE                                                           GL0892
      DATA NGAUSS/6/
C
C    FIND START AND END OF GAUSSIAN
C
      IFA=NZTYPE(NI)*4-3
      IF(C(IFA+1,1).NE.0.D0)THEN
         ILA=IFA+3
      ELSE
         ILA=IFA
      ENDIF
      IFB=NZTYPE(NJ)*4-3
      IF(C(IFB+1,1).NE.0.D0)THEN
         ILB=IFB+3
      ELSE
         ILB=IFB
      ENDIF
C
C  CONVERT R INTO AU
C
      R=R/0.529167D0
      R = R**2
      KA=0
      DO 80 I=IFA,ILA
         KA=KA+1
         NAT=KA-1
         KB=0
         DO 80 J=IFB,ILB
            KB=KB+1
            NBT=KB-1
C
C         DECIDE IS IT AN S-S, S-P, P-S, OR P-P OVERLAP
C
            IF(NAT.GT.0.AND.NBT.GT.0) THEN
C    P-P
               IS=4
               TOMB=(XI(NAT)-XJ(NAT))*(XI(NBT)
     1-XJ(NBT))*3.5711928576D0
            ELSEIF(NAT.GT.0) THEN
C    P-S
               IS=3
               TOMB=(XI(NAT)-XJ(NAT))*1.88976D0
            ELSEIF(NBT.GT.0) THEN
C    S-P
               IS=2
               TOMB=(XI(NBT)-XJ(NBT))*1.88976D0
            ELSE
C    S-S
               IS=1
            ENDIF
            DO 60 K=1,NGAUSS
               DO 60 L=1,NGAUSS
                  S(K,L)=0.0D0
                  AMB=Z(I,K)+Z(J,L)
                  APB=Z(I,K)*Z(J,L)
                  ADB=APB/AMB
C
C           CHECK OF OVERLAP IS NON-ZERO BEFORE STARTING
C
                  IF((ADB*R).LT.90.D0) THEN
                     ABN=1.0D0
                     GO TO(50,10,20,30),IS
   10                ABN=2.*TOMB*Z(I,K)*SQRT(Z(J,L))/AMB
                     GO TO 50
   20                ABN=-2.*TOMB*Z(J,L)*SQRT(Z(I,K))/AMB
                     GO TO 50
   30                ABN=-ADB*TOMB
                     IF(NAT.EQ.NBT) ABN=ABN+0.5D0
   40                ABN=4.0D0*ABN*SQRT(APB)/AMB
   50                S(K,L)=SQRT((2.*SQRT(APB)/AMB)**3)*EXP(-ADB*R)*ABN
                  ENDIF
   60       CONTINUE
            SG(KA,KB)=0.0D0
            DO 70 K=1,NGAUSS
               DO 70 L=1,NGAUSS
   70       SG(KA,KB)=SG(KA,KB)+S(K,L)*C(I,K)*C(J,L)
   80 CONTINUE
      RETURN
      END
************************************************************************
*
*  GRID CALCULATES THE ENERGY-SURFACE RESULTING FROM VARIATION OF
*       TWO COORDINATES. THE STEP-SIZE IS STEP1 AND STEP2, AND A 11
*       BY 11 GRID OF POINTS IS GENERATED
*
************************************************************************
      SUBROUTINE GRID
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /GEOVAR/ XPARAM(MAXPAR), NVAR, LOC(2,MAXPAR)               IR0394
      COMMON /GRADNT/ GRAD(MAXPAR),GNORM
      COMMON /GRAVEC/ COSINE
      COMMON /MESH  / LATOM1, LPARA1, LATOM2, LPARA2
      COMMON /KEYWRD/ KEYWRD
      DIMENSION SURFAC(11,11)
      CHARACTER*80 KEYWRD
         SAVE                                                           GL0892
      STEP1=READA(KEYWRD,INDEX(KEYWRD,'STEP1')+6)
      STEP2=READA(KEYWRD,INDEX(KEYWRD,'STEP2')+6)
C
C  THE CENTRAL VALUE OF THE FIRST AND SECOND DIMENSIONS ARE
C      GEO(LPARA1,LATOM1) AND GEO(LPARA2,LATOM2)
      NPTS=11
C NPTS MUST BE ODD, IN ORDER TO HAVE A CENTER POINT.
      NPTS2=NPTS/2
      DEGREE=180.D0/3.14159265359D0
      IF(LPARA1.NE.1)STEP1=STEP1/DEGREE
      IF(LPARA2.NE.1)STEP2=STEP2/DEGREE
      START1=GEO(LPARA1,LATOM1)-(NPTS2+1)*STEP1
      START2=GEO(LPARA2,LATOM2)-(NPTS2+1)*STEP2
C
C  NOW TO SWEEP THROUGH THE GRID OF POINTS LEFT TO RIGHT THEN RIGHT
C  TO LEFT. THIS SHOULD AVOID THE GEOMETRY OR SCF GETTING MESSED UP.
C
      GEO(LPARA1,LATOM1)=START1
      GEO(LPARA2,LATOM2)=START2
      IONE=-1.D0
      IF(LPARA1.NE.1) THEN
         C1=DEGREE
      ELSE
         C1=1.D0
      ENDIF
      IF(LPARA2.NE.1) THEN
         C2=DEGREE
      ELSE
         C2=1.D0
      ENDIF
      WRITE(6,'(''   FIRST VARIABLE   SECOND VARIABLE   FUNCTION'')')
      DO 20 ILOOP=1,NPTS
         GEO(LPARA1,LATOM1)=GEO(LPARA1,LATOM1)+STEP1
         IONE=-IONE
         JLOOP1=0
         IF(IONE.LT.0)JLOOP1=NPTS+1
         DO 10 JLOOP=1,NPTS
            JLOOP1=JLOOP1+IONE
            GEO(LPARA2,LATOM2)=GEO(LPARA2,LATOM2)+STEP2*IONE
            CALL FLEPO(XPARAM, NVAR, ESCF)
            SURFAC(ILOOP,JLOOP1)=ESCF
            WRITE(6,'('' :'',F16.5,F16.5,F13.6)')GEO(LPARA1,LATOM1)*C1,
     1        GEO(LPARA2,LATOM2)*C2,ESCF
   10    CONTINUE
         GEO(LPARA2,LATOM2)=GEO(LPARA2,LATOM2)+STEP2*IONE
   20 CONTINUE
      WRITE(6,'(/10X,''HORIZONTAL: VARYING SECOND PARAMETER,'',
     1          /10X,''VERTICAL:   VARYING FIRST PARAMETER'')')
      WRITE(6,'(/10X,''WHOLE OF GRID, SUITABLE FOR PLOTTING'',//)')
      DO 30 I=1,NPTS
   30 WRITE(6,'(11F7.2)')(SURFAC(J,I),J=1,NPTS)
      END
      SUBROUTINE HADDON (W,L,M,LOC,A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(3,*)
         SAVE                                                           GL0892
C**********************************************************************
C
C   HADDON CALCULATES THE VALUE OF A SYMMETRY-DEPENDENT VARIABLE
C
C  ON INPUT: M   = NUMBER SPECIFYING THE SYMMETRY OPERATION
C            LOC = ADDRESS OF REFERENCE ATOM
C            A   = ARRAY OF INTERNAL COORDINATES
C  ON OUTPUT W   = VALUE OF DEPENDENT FUNCTION
C            L   = 1 (FOR BOND LENGTH), 2 (ANGLE), OR 3 (DIHEDRAL)
C**********************************************************************
      PI = 3.1415926536D00
      IF (M.GT.18 .OR. M.LT.1) THEN
         WRITE(6,'(///10X,''UNDEFINED SYMMETRY FUNCTION USED'')')
         STOP
      ENDIF
      I=LOC
      GO TO
     1(140,160,10,20,30,40,50,60,70,80,90,100,110,120,150,170,180,190),
     2M
   10 W=A(3,I)
      GO TO 130
   20 W=(PI/2.0D00)-A(3,I)
      GO TO 130
   30 W=(PI/2.0D00)+A(3,I)
      GO TO 130
   40 W=(2.0D00*PI/3.0D00)-A(3,I)
      GO TO 130
   50 W=(2.0D00*PI/3.0D00)+A(3,I)
      GO TO 130
   60 W=(PI)-A(3,I)
      GO TO 130
   70 W=(PI)+A(3,I)
      GO TO 130
   80 W=(4.0D00*PI/3.0D00)-A(3,I)
      GO TO 130
   90 W=(4.0D00*PI/3.0D00)+A(3,I)
      GO TO 130
  100 W=(3.0D00*PI/2.0D00)-A(3,I)
      GO TO 130
  110 W=(3.0D00*PI/2.0D00)+A(3,I)
      GO TO 130
  120 W=-A(3,I)
  130 L=3
      RETURN
  140 L=1
      W=A(1,I)
      RETURN
  150 L=1
      W=A(1,I)/2.0D00
      RETURN
  160 L=2
      W=A(2,I)
      RETURN
  170 L=2
      W=A(2,I)/2.0D00
      RETURN
  180 L=2
      W=PI-A(2,I)
      RETURN
  190 CALL DEPVAR (A,I,W,L)
      RETURN
C
      END
      FUNCTION HELECT(N,P,H,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION P(*), H(*), F(*)
         SAVE                                                           GL0892
C***********************************************************************
C
C    SUBROUTINE CALCULATES THE ELECTRONIC ENERGY OF THE SYSTEM IN EV.
C
C    ON ENTRY N = NUMBER OF ATOMIC ORBITALS.
C             P = DENSITY MATRIX, PACKED, LOWER TRIANGLE.
C             H = ONE-ELECTRON MATRIX, PACKED, LOWER TRIANGLE.
C             F = TWO-ELECTRON MATRIX, PACKED, LOWER TRIANGLE.
C    ON EXIT
C        HELECT = ELECTRONIC ENERGY.
C
C    NO ARGUMENTS ARE CHANGED.
C
C***********************************************************************
      ED=0.0D00
      EE=0.0D00
      K=0
      NN=N+1
      DO 20 I=2,NN
         K=K+1
         JJ=I-1
         ED=ED+P(K)*(H(K)+F(K))
         IF (I.EQ.NN) GO TO 20
         DO 10 J=1,JJ
            K=K+1
   10    EE=EE+P(K)*(H(K)+F(K))
   20 CONTINUE
      EE=EE+.5D00*ED
      HELECT=EE
      RETURN
C
      END
      SUBROUTINE HQRII(A,N,M,E,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
************************************************************************
      DIMENSION A(*), E(N), V(N,M)
*************************************************************
*
* HQRII IS A DIAGONALISATION ROUTINE, WRITTEN BY YOSHITAKA BEPPU OF
*       NAGOYA UNIVERSITY, JAPAN.
*       FOR DETAILS SEE 'COMPUTERS & CHEMISTRY' VOL.6 1982. PAGE 000.
*
* ON INPUT    A       = MATRIX TO BE DIAGONALIZED
*             N       = SIZE OF MATRIX TO BE DIAGONALIZED.
*             M       = NUMBER OF EIGENVECTORS NEEDED.
*             E       = ARRAY OF SIZE AT LEAST N
*             V       = ARRAY OF SIZE AT LEAST NMAX*M
*
* ON OUTPUT   E       = EIGENVALUES
*             V       = EIGENVECTORS IN ARRAY OF SIZE NMAX*M
*
************************************************************************
      parameter(isize=3*numatm)
      dimension fv1(isize),sqa(isize*isize)
         SAVE                                                           GL0892
      l=1
      do 5 j=1,n
      do 5 i=1,j
      sqa(j+(i-1)*n)=a(l)
      l=l+1
5     continue
      call tred3(n,n,sqa,e,fv1,v)
      call tql2(n,n,e,fv1,v,ierr)
      RETURN
      END
      SUBROUTINE IJKL(I1, I2, J1, J2, ELEM, A1, MDIM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION A1(MDIM,NMECI)
      COMMON /XYIJKL/ XY(NMECI,NMECI,NMECI,NMECI)
************************************************************************
*
*   IJKL FILLS THE TWO-ELECTRON MATRIX XY WITH REPULSION INTEGRALS.
*        XY(I,J,K,L) IS THE REPULSION BETWEEN ONE ELECTRON IN
*        M.O.S I AND J AND AN ELECTRON IN M.O.S K AND L.
*        <I1(1),J1(1)/I2(2),J2(2)>
************************************************************************
      COMMON /WMATRX/ WJ(N2ELEC), WK(N2ELEC)
      DIMENSION W(N2ELEC)
      EQUIVALENCE (W,WJ)
         SAVE                                                           GL0892
      IF (XY(I1,J1,I2,J2).EQ.100.0D0) THEN
         X= SPCG(A1(1,I1),A1(1,J1),A1(1,I2),A1(1,J2),W,WJ)
C#          WRITE(6,'(4I6,F13.6)')I1,J1,I2,J2,X
         XY(I1,J1,I2,J2)=X
         XY(I1,J1,J2,I2)=X
         XY(J1,I1,I2,J2)=X
         XY(J1,I1,J2,I2)=X
         XY(I2,J2,I1,J1)=X
         XY(I2,J2,J1,I1)=X
         XY(J2,I2,I1,J1)=X
         XY(J2,I2,J1,I1)=X
      ENDIF
      IF (XY(I1,I2,J1,J2).EQ.100.0D0) THEN
         Z= SPCG(A1(1,I1),A1(1,I2),A1(1,J1),A1(1,J2),W,WJ)
C#          WRITE(6,'(4I6,F13.6)')I1,I2,J1,J2,Z
         XY(I1,I2,J1,J2)=Z
         XY(I1,I2,J2,J1)=Z
         XY(I2,I1,J1,J2)=Z
         XY(I2,I1,J2,J1)=Z
         XY(J1,J2,I1,I2)=Z
         XY(J1,J2,I2,I1)=Z
         XY(J2,J1,I1,I2)=Z
         XY(J2,J1,I2,I1)=Z
      ENDIF
      IF (XY(I1,J2,I2,J1).EQ.100.0D0) THEN
         Y= SPCG(A1(1,I1),A1(1,J2),A1(1,I2),A1(1,J1),W,WJ)
C#          WRITE(6,'(4I6,F13.6)')I1,J2,I2,J1,Y
         XY(I1,J2,I2,J1)=Y
         XY(I1,J2,J1,I2)=Y
         XY(J2,I1,I2,J1)=Y
         XY(J2,I1,J1,I2)=Y
         XY(I2,J1,I1,J2)=Y
         XY(I2,J1,J2,I1)=Y
         XY(J1,I2,I1,J2)=Y
         XY(J1,I2,J2,I1)=Y
      ENDIF
      X=XY(I1,J1,I2,J2)
      Y=XY(I1,J2,J1,I2)
      ELEM=X-Y
      RETURN
      END
      SUBROUTINE ITER  (H, W, WJ, WK, EE, FULSCF,RAND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DOUBLE PRECISION MECI
cmgc
      common /cory/ e_cory
cmgc
      DIMENSION H(*), W(*), WJ(*), WK(*)
      COMMON /FOKMAT/ F(MPACK), FB(MPACK)
      COMMON /DENSTY/ P(MPACK), PA(MPACK), PB(MPACK)
      COMMON /VECTOR/ C(MORB2),EIGS(MAXORB),CBETA(MORB2),EIGB(MAXORB)
      COMMON /GRADNT/ DUMY(MAXPAR),GNORM
      COMMON /LAST  / LAST
      COMMON /MESAGE/ IFLEPO,IITER
      COMMON /ATHEAT/ ATHEAT
      COMMON /ENUCLR/ ENUCLR
      COMMON /CITERM/ XI,XJ,XK
C    changed common path for portability  (IR)
      COMMON /PATHI / LATOM,LPARAM
      COMMON /PATHR / REACT(200)
c
      COMMON /NUMCAL/ NUMCAL
      COMMON /TIMER / TIME0
      LOGICAL FULSCF, RAND
C***********************************************************************
C
C     ITER GENERATES A SCF FIELD AND RETURNS THE ENERGY IN "ENERGY"
C
C THE MAIN ARRAYS USED IN ITER ARE:
C            P      ONLY EVER CONTAINS THE TOTAL DENSITY MATRIX
C            PA     ONLY EVER CONTAINS THE ALPHA DENSITY MATRIX
C            PB     ONLY EVER CONTAINS THE BETA DENSITY MATRIX
C            C      ONLY EVER CONTAINS THE EIGENVECTORS
C            H      ONLY EVER CONTAINS THE ONE-ELECTRON MATRIX
C            F      STARTS OFF CONTAINING THE ONE-ELECTRON MATRIX,
C                   AND IS USED TO HOLD THE FOCK MATRIX
C            W      ONLY EVER CONTAINS THE TWO-ELECTRON MATRIX
C
C THE MAIN INTEGERS CONSTANTS IN ITER ARE:
C
C            LINEAR SIZE OF PACKED TRIANGLE = NORBS*(NORBS+1)/2
C
C THE MAIN INTEGER VARIABLES ARE
C            NITER  NUMBER OF ITERATIONS EXECUTED
C
C  PRINCIPAL REFERENCES:
C
C   ON MNDO: "GROUND STATES OF MOLECULES. 38. THE MNDO METHOD.
C             APPROXIMATIONS AND PARAMETERS."
C             DEWAR, M.J.S., THIEL,W., J. AM. CHEM. SOC.,99,4899,(1977).
C   ON SHIFT: "THE DYNAMIC 'LEVEL SHIFT' METHOD FOR IMPROVING THE
C             CONVERGENCE OF THE SCF PROCEDURE", A. V. MITIN, J. COMP.
C             CHEM. 9, 107-110 (1988)
C   ON HALF-ELECTRON: "MINDO/3 COMPARISON OF THE GENERALIZED S.C.F.
C             COUPLING OPERATOR AND "HALF-ELECTRON" METHODS FOR
C             CALCULATING THE ENERGIES AND GEOMETRIES OF OPEN SHELL
C             SYSTEMS"
C             DEWAR, M.J.S., OLIVELLA, S., J. CHEM. SOC. FARA. II,
C             75,829,(1979).
C   ON PULAY'S CONVERGER: "CONVERGANCE ACCELERATION OF ITERATIVE
C             SEQUENCES. THE CASE OF SCF ITERATION", PULAY, P.,
C             CHEM. PHYS. LETT, 73, 393, (1980).
C   ON CNVG:  IT ENCORPORATES THE IMPROVED ITERATION SCHEME (IIS) BY
C             PIOTR BADZIAG & FRITZ SOLMS. ACCEPTED FOR PUBLISHING
C             IN COMPUTERS & CHEMISTRY
C   ON PSEUDODIAGONALISATION: "FAST SEMIEMPIRICAL CALCULATIONS",
C             STEWART. J.J.P., CSASZAR, P., PULAY, P., J. COMP. CHEM.,
C             3, 227, (1982)
C
C***********************************************************************
      DIMENSION POLD(MPACK), POLD2(MPACK), POLD3(MAXORB+400)
      DIMENSION PBOLD(MPACK), PBOLD2(MPACK), PBOLD3(MAXORB+400)
************************************************************************
*                                                                      *
*   IF, FOR ANY REASON, ITER FAILS TO GENERATE AN SCF, THEN UNCOMMENT  *
*   ALL LINES THAT BEGIN WITH A ''.  THE EFFECT OF DOING THIS IS TO  *
*   INCREASE THE STORAGE REQUIREMENT BY ABOUT 20% AND TO ALLOW PULAY'S *
*   AND CAMP AND KING'S CONVERGERS TO BE USED.                         *
*                                                                      *
************************************************************************
      DIMENSION  AR1(2*NPULAY), AR2(2*NPULAY), AR3(2*NPULAY),
     1 AR4(2*NPULAY)
      DIMENSION  BR1(2*NPULAY), BR2(2*NPULAY), BR3(2*NPULAY),
     1 BR4(2*NPULAY)
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /MOLORB/ DUMMY(MAXORB),PDIAG(MAXORB)
      COMMON /KEYWRD/ KEYWRD
      COMMON /NUMSCF/ NSCF
      CHARACTER KEYWRD*80, ABPRT(3)*5
      LOGICAL PRTFOK,PRTEIG,PRTDEN, DEBUG, PRTENG, TIMES, CI
     1,UHF, NEWDG, SCF1, HALFE, FORCE, PRT1EL,PRTPL, OKNEWD
     2,EXCITD, MINPRT, FRST, BFRST, OKPULY, READY, PRTVEC,
     3CAMKIN, ALLCON, MAKEA, MAKEB, INCITR, CAPPS
         SAVE                                                           GL0892
      DATA ICALCN/0/, DEBUG/.FALSE./, PRTFOK/.FALSE./
      DATA PRTEIG/.FALSE./,PRTDEN,PRTENG/.FALSE.,.FALSE./
      DATA PRT1EL/.FALSE./
      DATA ABPRT/'     ','ALPHA',' BETA'/
C
C  INITIALIZE
C
      IFIL=0
      IHOMO=NCLOSE+NALPHA
      IHOMOB=NCLOSE+NBETA
      EOLD=1.D2
      READY=.FALSE.
      IF (ICALCN.NE.NUMCAL) THEN
         SHIFT=0.D0
         ICALCN=NUMCAL
         SHFMAX=20.D0
         LINEAR=(NORBS*(NORBS+1))/2
C
C    DEBUG KEY-WORDS WORKED OUT
C
         MINPRT=(INDEX(KEYWRD,'SADDLE')+
     1      LATOM .EQ.0 .OR. DEBUG)
         PRTEIG=( INDEX(KEYWRD,'EIGS') .NE. 0 )
         PRTENG=( INDEX(KEYWRD,'ENERGY').NE.0 )
         PRTPL =( INDEX(KEYWRD,' PL ')  .NE.0 )
         DEBUG=( INDEX(KEYWRD,'DEBUG') .NE. 0 )
         PRT1EL=( INDEX(KEYWRD,'1ELEC') .NE.0 .AND. DEBUG)
         PRTDEN=( INDEX(KEYWRD,'DENS').NE.0 .AND. DEBUG)
         PRTFOK=( INDEX(KEYWRD,'FOCK') .NE. 0  .AND. DEBUG)
         PRTVEC=( INDEX(KEYWRD,'VECT') .NE. 0  .AND. DEBUG)
         DEBUG=( INDEX(KEYWRD,'ITER') .NE. 0 )
C
C INITIALIZE SOME LOGICALS AND CONSTANTS
C
         NEWDG =.FALSE.
         PL    =1.D0
      bshift = 0.0
         SHIFT=1.D0
*
* SCFCRT AND PLTEST ARE MACHINE-PRECISION DEPENDENT
*
         SCFCRT=1.D-6
         PLTEST=0.0001D0
         ITRMAX = 200
         NMOS=0
         NCIS=0
         NA2EL=NCLOSE
         NA1EL=NALPHA+NOPEN
         NB2EL=0
         NB1EL=NBETA+NOPEN
C
C  USE KEY-WORDS TO ASSIGN VARIOUS CONSTANTS
C
         IF(INDEX(KEYWRD,'C.I.').NE.0)
     1      NMOS=READA(KEYWRD,INDEX(KEYWRD,'C.I.')+5)
         IF(INDEX(KEYWRD,'MICROS').NE.0)
     1      NCIS=READA(KEYWRD,INDEX(KEYWRD,'MICROS'))
         IF(INDEX(KEYWRD,'FILL').NE.0)
     1      IFILL=-READA(KEYWRD,INDEX(KEYWRD,'FILL'))
         IF(INDEX(KEYWRD,'SHIFT').NE.0)
     1      BSHIFT=-READA(KEYWRD,INDEX(KEYWRD,'SHIFT'))
         IF(BSHIFT.NE.0)TEN=BSHIFT
         IF(INDEX(KEYWRD,'ITRY').NE.0)
     1      ITRMAX=READA(KEYWRD,INDEX(KEYWRD,'ITRY'))
         CAMKIN=(INDEX(KEYWRD,'KING')+INDEX(KEYWRD,'CAMP') .NE. 0)
         CI=(INDEX(KEYWRD,'MICROS')+INDEX(KEYWRD,'C.I.') .NE. 0)
         OKPULY=.FALSE.
         OKPULY=(INDEX(KEYWRD,'PULAY').NE.0)
         UHF=(INDEX(KEYWRD,'UHF') .NE. 0)
         SCF1=(INDEX(KEYWRD,'1SCF') .NE. 0)
         OKNEWD=ABS(BSHIFT) .LT. 0.001D0
         IF(CAMKIN.AND.ABS(BSHIFT).GT.1.D-5) BSHIFT=4.44D0
         EXCITD=(INDEX(KEYWRD,'EXCITED') .NE. 0)
         TIMES=(INDEX(KEYWRD,'TIMES') .NE. 0)
         FORCE=(INDEX(KEYWRD,'FORCE') .NE. 0)
         ALLCON=(OKPULY.OR.CAMKIN)
C
C   SET UP C.I. PARAMETERS
C   NMOS IS NO. OF M.O.S USED IN C.I.
C   NCIS IS CHANGE IN SPIN, OR NUMBER OF STATES
C
         IF(NMOS.EQ.0) NMOS=NOPEN-NCLOSE
         IF(NCIS.EQ.0) THEN
            IF(INDEX(KEYWRD,'TRIPLET')+INDEX(KEYWRD,'QUART').NE.0)NCIS=1
            IF(INDEX(KEYWRD,'QUINTET')+INDEX(KEYWRD,'SEXTE').NE.0)NCIS=2
         ENDIF
C
C   DO WE NEED A CAPPED ATOM CORRECTION?
C
         CAPPS=.FALSE.
         DO 10 I=1,NUMAT
   10    IF(NAT(I).EQ.102)CAPPS=.TRUE.
         IITER=1
         TRANS=0.1D0
         IF(INDEX(KEYWRD,'RESTART')+INDEX(KEYWRD,'OLDENS')
     1      .NE. 0) THEN
            IF(INDEX(KEYWRD,'OLDENS').NE.0)
     1   OPEN(UNIT=10,FILE='FOR010',STATUS='UNKNOWN',FORM='UNFORMATTED')
            REWIND 10
            READ(10)(PA(I),I=1,LINEAR)
            IF( UHF) THEN
               READ(10)(PB(I),I=1,LINEAR)
               DO 20 I=1,LINEAR
                  POLD(I)=PA(I)
                  PBOLD(I)=PB(I)
   20          P(I)=PA(I)+PB(I)
            ELSE
               DO 30 I=1,LINEAR
                  PB(I)=PA(I)
                  PBOLD(I)=PA(I)
                  POLD(I)=PA(I)
   30          P(I)=PA(I)*2.D0
            ENDIF
         ELSE
            NSCF=0
            DO 40 I=1,LINEAR
               P(I)=0.D0
               PA(I)=0.D0
   40       PB(I)=0.D0
            W1=NA1EL/(NA1EL+1.D-6+NB1EL)
            W2=1.D0-W1
            IF(W1.LT.1.D-6)W1=0.5D0
            IF(W2.LT.1.D-6)W2=0.5D0
            RANDOM=1.0D0
            IF(UHF.AND.NA1EL.EQ.NB1EL) RANDOM=1.1D0
            IF(.NOT.RAND)RANDOM=0.D0
            DO 50 I=1,NORBS
               J=(I*(I+1))/2
               P(J)=PDIAG(I)
               PA(J)=P(J)*W1*RANDOM
               RANDOM=1.D0/RANDOM
   50       PB(J)=P(J)*W2*RANDOM
            DO 60 I=1,LINEAR
               PBOLD(I)=PB(I)
   60       POLD(I)=PA(I)
         ENDIF
         HALFE=(NOPEN .NE. NCLOSE)
C
C   DETERMINE THE SELF-CONSISTENCY CRITERION
c---Note:  scfcrt is a machine dependent constant
c---set to 1.0d-12 for vax machines
c---set to 1.0d-10 for Cray machines
C
         IF( HALFE .OR.  INDEX(KEYWRD,'PREC') .NE. 0)
     3                               SCFCRT=SCFCRT*0.01D0
         IF( INDEX(KEYWRD,'POLAR') + INDEX(KEYWRD,'NLLSQ') +
     2 INDEX(KEYWRD,'SIGMA') .NE. 0) SCFCRT=SCFCRT*0.01D0
         IF(FORCE)                   SCFCRT=SCFCRT*0.01D0
         SCFCRT=MAX(SCFCRT,1.D-10)
C
C  THE USER CAN STATE THE SCF CRITERION, IF DESIRED.
C
         I=INDEX(KEYWRD,'SCFCRT')
         IF(I.NE.0) THEN
            SCFCRT=READA(KEYWRD,I)
            WRITE(6,'(''  SCF CRITERION ='',G14.4)')SCFCRT
            IF(SCFCRT.LT.1.D-10)
     1 WRITE(6,'(//2X,'' THERE IS A RISK OF INFINITE LOOPING WITH'',
     2'' THE SCFCRT LESS THAN 1.D-10'')')
         ELSE
            IF(DEBUG)WRITE(6,'(''  SCF CRITERION ='',G13.9)')SCFCRT
         ENDIF
         LAST=0
C
C   END OF INITIALIZATION SECTION.
C
      ELSEIF(FORCE)THEN
C
C   RESET THE DENSITY MATRIX IF MECI HAS FORMED AN EXCITED STATE
C   SUM IS NOT USED
C
         SUM=MECI(EIGS,C,CBETA,EIGB, NORBS,NMOS,1,.FALSE.)
      ENDIF
C
C   INITIALIZATION OPERATIONS DONE EVERY TIME ITER IS CALLED
C
      MAKEA=.TRUE.
      MAKEB=.TRUE.
      IF(NEWDG) NEWDG=(ABS(BSHIFT).LT.0.001D0)
      IF(LAST.EQ.1) NEWDG=.FALSE.
      SELCON=SCFCRT*23.061
      IF(PLTEST.LT.SCFCRT) PLTEST=SCFCRT
      IF(NALPHA.NE.NBETA.OR..NOT.UHF)PLTEST=0.001D0
      IF(.NOT. FORCE .AND. .NOT. HALFE) THEN
         IF(GNORM.GT.5.D0) SELCON=SCFCRT*GNORM*0.2D0
         IF(GNORM.GT.200.D0) SELCON=SCFCRT*50.D0
      ENDIF
      IF(DEBUG)WRITE(6,'(''  SELCON, GNORM'',2F16.7)')SELCON,GNORM
cmgc      TITER1=SECOND()
      IF(PRT1EL) THEN
         WRITE(6,'(//10X,''ONE-ELECTRON MATRIX AT ENTRANCE TO ITER'')')
         CALL VECPRT(H,NORBS)
      ENDIF
      IREDY=1
   70 NITER=0
cmgc      TIME1=SECOND()
      FRST=.TRUE.
      IF(CAMKIN) THEN
         MODEA=1
         MODEB=1
      ELSE
         MODEA=0
         MODEB=0
      ENDIF
      BFRST=.TRUE.
**********************************************************************
*                                                                    *
*                                                                    *
*                START THE SCF LOOP HERE                             *
*                                                                    *
*                                                                    *
**********************************************************************
      INCITR=.TRUE.
   80 INCITR=(MODEA.NE.3.AND.MODEB.NE.3)
      IF(INCITR)NITER=NITER+1
      IF(NITER.GT.ITRMAX-10.AND..NOT.ALLCON) THEN
************************************************************************
*                                                                      *
*                   SWITCH ON ALL CONVERGERS                           *
*                                                                      *
************************************************************************
         WRITE(6,'(//,'' ALL CONVERGERS ARE NOW FORCED ON'',/
     1          '' SHIFT=10, PULAY ON, CAMP-KING ON'',/
     2          '' AND ITERATION COUNTER RESET'',//)')
         ALLCON=.TRUE.
         BSHIFT=4.44D0
         IREDY=-4
         EOLD=100.D0
         OKPULY=.TRUE.
         NEWDG=.FALSE.
         CAMKIN=(.NOT.HALFE)
         GOTO 70
      ENDIF
************************************************************************
*                                                                      *
*                        MAKE THE ALPHA FOCK MATRIX                    *
*                                                                      *
************************************************************************
      IF(ABS(SHIFT).GT.1.D-10.AND.BSHIFT .NE. 0.D0) THEN
         L=0
         IF(NITER.GT.1)THEN
            IF(NEWDG.AND..NOT.(HALFE.OR.CAMKIN))THEN
C
C  SHIFT WILL APPLY TO THE VIRTUAL ENERGY LEVELS USED IN THE
C  PSEUDODIAGONALIIZATION. IF DIFF IS -VE, GOOD, THEN LOWER THE
C  HOMO-LUMO GAP BY 0.1EV, OTHERWISE INCREASE IT.
               IF(DIFF.GT.0)THEN
               SHIFT=1.D0
C
C IF THE PSEUDODIAGONALIZATION APPROXIMATION -- THAT THE WAVEFUNCTION
C IS ALMOST STABLE -- IS INVALID, TURN OFF NEWDG
               IF(DIFF.GT.1)NEWDG=.FALSE.
               ELSE
               SHIFT=-0.1D0
               ENDIF
            ELSE
               SHIFT=TEN+EIGS(IHOMO+1)-EIGS(IHOMO)+SHIFT
            ENDIF
            IF(DIFF.GT.0.D0) THEN
               IF(SHIFT.GT.4.D0)SHFMAX=4.5D0
               IF(SHIFT.GT.SHFMAX)SHFMAX=MAX(SHFMAX-0.5D0,0.D0)
            ENDIF
C
C   IF SYSTEM GOES UNSTABLE, LIMIT SHIFT TO THE RANGE -INFINITY - SHFMAX
C   BUT IF SYSTEM IS STABLE, LIMIT SHIFT TO THE RANGE -INFINITY - +20
C
            SHIFT=MAX(-20.D0,MIN(SHFMAX,SHIFT))
            IF(ABS(SHIFT-SHFMAX).LT.1.D-5)SHFMAX=SHFMAX+0.01D0
            IF(ABS(BSHIFT-4.44D0).LT.1.D-5)SHIFT=-8.D0
            IF(UHF)THEN
               IF(NEWDG.AND..NOT.(HALFE.OR.CAMKIN))THEN
                  SHIFTB=TEN-TENOLD
               ELSE
                  SHIFTB=TEN+EIGS(IHOMOB+1)-EIGS(IHOMOB)+SHIFTB
               ENDIF
               IF(DIFF.GT.0.D0)SHIFTB=MIN(4.D0,SHIFTB)
               SHIFTB=MAX(-20.D0,MIN(SHFMAX,SHIFTB))
               IF(ABS(BSHIFT-4.44D0).LT.1.D-5)SHIFTB=-8.D0
C#       WRITE(6,*)'SHIFT:',SHIFT,SHIFTB
               DO 90 I=IHOMOB+1,NORBS
   90          EIGB(I)=EIGB(I)+SHIFTB
            ELSE
C#       WRITE(6,*)'SHIFT:',SHIFT
            ENDIF
         ENDIF
         TENOLD=TEN
         DO 100 I=IHOMO+1,NORBS
  100    EIGS(I)=EIGS(I)+SHIFT
         DO 120 I=1,NORBS
            DO 110 J=1,I
               L=L+1
  110       F(L)=H(L)+SHIFT*PA(L)
  120    F(L)=F(L)-SHIFT
      ELSEIF (RAND.AND.LAST.EQ.0.AND.NITER.LT.2.AND.FULSCF)THEN
         RANDOM=0.001D0
         DO 130 I=1,LINEAR
            RANDOM=-RANDOM
  130    F(I)=H(I)+RANDOM
      ELSE
         DO 140 I=1,LINEAR
  140    F(I)=H(I)
      ENDIF
  150 CALL FOCK2(F,P,PA,W, WJ, WK,NUMAT,NFIRST,NMIDLE,NLAST)
      CALL FOCK1(F,P,PA,PB)
************************************************************************
*                                                                      *
*                        MAKE THE BETA FOCK MATRIX                     *
*                                                                      *
************************************************************************
      IF (UHF) THEN
         IF(SHIFTB .NE. 0.D0) THEN
            L=0
            DO 170 I=1,NORBS
               DO 160 J=1,I
                  L=L+1
  160          FB(L)=H(L)+SHIFTB*PB(L)
  170       FB(L)=FB(L)-SHIFTB
         ELSEIF (RAND.AND.LAST.EQ.0.AND.NITER.LT.2.AND.FULSCF)THEN
            RANDOM=0.001
            DO 180 I=1,LINEAR
               RANDOM=-RANDOM
  180       FB(I)=H(I)+RANDOM
         ELSE
            DO 190 I=1,LINEAR
  190       FB(I)=H(I)
         ENDIF
         CALL FOCK2(FB,P,PB,W, WJ, WK,NUMAT,NFIRST,NMIDLE,NLAST)
         CALL FOCK1(FB,P,PB,PA)
      ENDIF
      IF( .NOT. FULSCF) GOTO 290
      IF(PRTFOK) THEN
         WRITE(6,200)NITER
  200    FORMAT('   FOCK MATRIX ON ITERATION',I3)
         CALL VECPRT (F,NORBS)
      ENDIF
************************************************************************
*                                                                      *
*                        CALCULATE THE ENERGY IN KCAL/MOLE             *
*                                                                      *
************************************************************************
      EE=HELECT(NORBS,PA,H,F)
      IF(UHF)THEN
         EE=EE+HELECT(NORBS,PB,H,FB)
      ELSE
         EE=EE*2.D0
      ENDIF
      IF(CAPPS)EE=EE+CAPCOR(NAT,NFIRST,NLAST,NUMAT,P,H)
      SCORR=SHIFT*(NOPEN-NCLOSE)*23.061D0*0.25D0*(FRACT*(2.D0-FRACT))
      ESCF=(EE+ENUCLR)*23.061D0+ATHEAT+SCORR
      IF(INCITR)THEN
         DIFF=ESCF-EOLD
         IF(DIFF.GT.0)THEN
            TEN=TEN-1.D0
         ELSE
            TEN=TEN*0.975D0+0.05D0
         ENDIF
C#         WRITE(6,'(2F12.6)')TEN,DIFF
         IF( (NITER.GT.4).AND. (PL.LT.PLTEST) .AND.                       IR0295
     &       (ABS(DIFF/MAX(1.D0,ABS(EE))).LT.SELCON) .AND. READY) THEN    
************************************************************************
*                                                                      *
*          SELF-CONSISTENCY TEST, EXIT MODE FROM ITERATIONS            *
*                                                                      *
************************************************************************
            IF (ABS(SHIFT) .LT. 1.D-10) GOTO 290
            SHIFT=0.D0
            SHIFTB=0.D0
            DO 210 I=1,LINEAR
  210       F(I)=H(I)
            MAKEA=.TRUE.
            MAKEB=.TRUE.
            GOTO 150
         ENDIF
        READY=(IREDY.GT.0.AND.(ABS(DIFF).LT.SELCON*10.D0.OR.PL.EQ.0.D0))
        IREDY=IREDY+1
      ENDIF
      IF(PRTPL.OR.DEBUG.AND.NITER.GT.ITRMAX-20) THEN
         IF(ABS(ESCF).GT.99999.D0) ESCF=SIGN(9999.D0,ESCF)
         IF(ABS(DIFF).GT.9999.D0)DIFF=0.D0
         IF(INCITR)
     1    WRITE(6,'('' ITERATION'',I3,'' PLS='',2E10.3,'' ENERGY  '',
     2F14.7,'' DELTAE'',F13.7)')NITER,PL,PLB,ESCF,DIFF
      ENDIF
      IF(INCITR)EOLD=ESCF
************************************************************************
*                                                                      *
*                        INVOKE THE CAMP-KING CONVERGER                *
*                                                                      *
************************************************************************
      IF(NITER.GT.2 .AND. CAMKIN .AND. MAKEA)
     1CALL INTERP(NORBS,NA1EL,NORBS-NA1EL, MODEA, ESCF/23.061D0,
     2F, C, AR1, AR2, AR3, AR4, AR1)
      MAKEB=.FALSE.
      IF(MODEA.EQ.3)GOTO 230
      MAKEB=.TRUE.
      IF( NEWDG ) THEN
************************************************************************
*                                                                      *
*                        INVOKE PULAY'S CONVERGER                      *
*                                                                      *
************************************************************************
         IF(OKPULY.AND.MAKEA.AND.IREDY.GT.1)
     1CALL PULAY(F,PA,NORBS,POLD,POLD2,POLD3,JALP,IALP,MPACK,FRST,PL)
************************************************************************
*                                                                      *
*           DIAGONALIZE THE ALPHA OR RHF SECULAR DETERMINANT           *
* WHERE POSSIBLE, USE THE PULAY-STEWART METHOD, OTHERWISE USE BEPPU'S  *
*                                                                      *
************************************************************************
         IF (HALFE.OR.CAMKIN) THEN
            CALL RSP(F,NORBS,NORBS,EIGS,C)
         ELSE
            CALL DIAG (F,C,NA1EL,EIGS,NORBS,NORBS)
         ENDIF
      ELSE
         CALL RSP(F,NORBS,NORBS,EIGS,C)
      ENDIF
      J=1
      IF(PRTVEC) THEN
         J=1
         IF(UHF)J=2
         WRITE(6,'(//10X,A,
     1'' EIGENVECTORS AND EIGENVALUES ON ITERATION'',I3)')
     2   ABPRT(J),NITER
         CALL MATOUT(C,EIGS,NORBS,NORBS,NORBS)
      ELSE
         IF (PRTEIG) WRITE(6,220)ABPRT(J),NITER,(EIGS(I),I=1,NORBS)
      ENDIF
  220 FORMAT(10X,A,'  EIGENVALUES ON ITERATION',I3,/10(6G13.6,/))
  230 IF(IFILL.NE.0)CALL SWAP(C,NORBS,NORBS,NA2EL,IFILL)
************************************************************************
*                                                                      *
*            CALCULATE THE ALPHA OR RHF DENSITY MATRIX                 *
*                                                                      *
************************************************************************
      IF(UHF)THEN
         CALL DENSIT( C,NORBS, NORBS, NA2EL,NA1EL, FRACT, PA, 1)
         IF(MODEA.NE.3.AND..NOT. (NEWDG.AND.OKPULY))
     1    CALL CNVG(PA, POLD, POLD2, NORBS, NITER, PL)
      ELSE
         CALL DENSIT( C,NORBS, NORBS, NA2EL,NA1EL, FRACT, P, 1)
         IF(MODEA.NE.3.AND..NOT. (NEWDG.AND.OKPULY))
     1    CALL CNVG(P, POLD, POLD2, NORBS, NITER, PL)
      ENDIF
************************************************************************
*                                                                      *
*                       UHF-SPECIFIC CODE                              *
*                                                                      *
************************************************************************
      IF( UHF )THEN
************************************************************************
*                                                                      *
*                        INVOKE THE CAMP-KING CONVERGER                *
*                                                                      *
************************************************************************
         IF(NITER.GT.2 .AND. CAMKIN .AND. MAKEB )
     1CALL INTERP(NORBS,NB1EL,NORBS-NB1EL, MODEB, ESCF/23.061D0,
     2FB, CBETA, BR1, BR2, BR3, BR4, BR1)
         MAKEA=.FALSE.
         IF(MODEB.EQ.3) GOTO 240
         MAKEA=.TRUE.
         IF( NEWDG ) THEN
************************************************************************
*                                                                      *
*                        INVOKE PULAY'S CONVERGER                      *
*                                                                      *
************************************************************************
            IF( OKPULY.AND.MAKEB.AND.IREDY.GT.1)
     1CALL PULAY(FB,PB,NORBS,PBOLD,PBOLD2,
     2PBOLD3,JBET,IBET,MPACK,BFRST,PLB)
************************************************************************
*                                                                      *
*           DIAGONALIZE THE ALPHA OR RHF SECULAR DETERMINANT           *
* WHERE POSSIBLE, USE THE PULAY-STEWART METHOD, OTHERWISE USE BEPPU'S  *
*                                                                      *
************************************************************************
            IF (HALFE.OR.CAMKIN) THEN
               CALL RSP(FB,NORBS,NORBS,EIGB,CBETA)
            ELSE
               CALL DIAG (FB,CBETA,NB1EL,EIGB,NORBS,NORBS)
            ENDIF
         ELSE
            CALL RSP(FB,NORBS,NORBS,EIGB,CBETA)
         ENDIF
         IF(PRTVEC) THEN
            WRITE(6,'(//10X,A,'' EIGENVECTORS AND EIGENVALUES ON '',
     1''ITERATION'',I3)')ABPRT(3),NITER
            CALL MATOUT(CBETA,EIGB,NORBS,NORBS,NORBS)
         ELSE
            IF (PRTEIG) WRITE(6,220)ABPRT(3),NITER,(EIGB(I),I=1,NORBS)
         ENDIF
************************************************************************
*                                                                      *
*                CALCULATE THE BETA DENSITY MATRIX                     *
*                                                                      *
************************************************************************
  240    CALL DENSIT( CBETA,NORBS, NORBS, NB2EL, NB1EL, FRACT, PB, 1)
         IF( .NOT. (NEWDG.AND.OKPULY))
     1CALL CNVG(PB, PBOLD, PBOLD2, NORBS, NITER, PLB)
      ENDIF
************************************************************************
*                                                                      *
*                   CALCULATE THE TOTAL DENSITY MATRIX                 *
*                                                                      *
************************************************************************
      IF(UHF) THEN
         DO 250 I=1,LINEAR
  250    P(I)=PA(I)+PB(I)
      ELSE
         DO 260 I=1,LINEAR
            PA(I)=P(I)*0.5D0
  260    PB(I)=PA(I)
      ENDIF
      IF(PRTDEN) THEN
         WRITE(6,'('' DENSITY MATRIX ON ITERATION'',I4)')NITER
         CALL VECPRT (P,NORBS)
      ENDIF
      OKNEWD=(PL.LT.SELCON .OR. OKNEWD)
      NEWDG=(PL.LT.TRANS .AND. OKNEWD .OR. NEWDG)
      IF(PL.LT.TRANS*0.3333D0)OKNEWD=.TRUE.
      IF (NITER .GE. ITRMAX) THEN
         IF(DIFF.LT.1.D-3.AND.PL.LT.1.D-4.AND..NOT.FORCE)THEN
            WRITE(6,'('' """""""""""""""UNABLE TO ACHIEVE SELF-CONSISTEN
     1CE, JOB CONTINUING'')')
            GOTO 290
         ENDIF
         IF(MINPRT)WRITE (6,270)
  270    FORMAT (//10X,'"""""""""""""UNABLE TO ACHIEVE SELF-CONSISTENCE'
     1,/)
         WRITE (6,280) DIFF,PL
  280    FORMAT (//,10X,'DELTAE= ',E12.4,5X,'DELTAP= ',E12.4,///)
         IFLEPO=9
         IITER=2
         CALL WRITMO (TIME0,ESCF)
C GA-FRIENDLY : change the next STOP to RETURN if run with GA code
C        RETURN                                                          GA-ON
         STOP
      ENDIF
      GO TO 80
**********************************************************************
*                                                                    *
*                                                                    *
*                      END THE SCF LOOP HERE                         *
*                NOW CALCULATE THE ELECTRONIC ENERGY                 *
*                                                                    *
*                                                                    *
**********************************************************************
*          SELF-CONSISTENCE ACHEIVED.
*
  290 EE=HELECT(NORBS,PA,H,F)
      IF(UHF) THEN
         EE=EE+HELECT(NORBS,PB,H,FB)
      ELSE
         EE=EE*2.D0 +
     1SHIFT*(NOPEN-NCLOSE)*23.061D0*0.25D0*(FRACT*(2.D0-FRACT))
      ENDIF
      IF(CAPPS)EE=EE+CAPCOR(NAT,NFIRST,NLAST,NUMAT,P,H)
C
C   NORMALLY THE EIGENVALUES ARE INCORRECT BECAUSE THE
C   PSEUDODIAGONALIZATION HAS BEEN USED.  IF THIS
C   IS THE LAST SCF, THEN DO AN EXACT DIAGONALIZATION
      IF(NSCF .EQ. 0 .OR. LAST .EQ. 1 .OR. CI .OR. HALFE) THEN
C
C  PUT F AND FB INTO POLD IN ORDER TO NOT DESTROY F AND FB
C  AND DO EXACT DIAGONALISATIONS
         DO 300 I=1,LINEAR
  300    POLD(I)=F(I)
         CALL RSP(POLD,NORBS,NORBS,EIGS,C)
         IF(UHF) THEN
            DO 310 I=1,LINEAR
  310       POLD(I)=FB(I)
            CALL RSP(POLD,NORBS,NORBS,EIGB,CBETA)
         ENDIF
         DO 320 I=1,LINEAR
  320    POLD(I)=P(I)
         IF(CI.OR.HALFE) THEN
            SUM=MECI(EIGS,C,CBETA,EIGB, NORBS,NMOS,0,.FALSE.)
            EE=EE+SUM
            IF(PRTPL)THEN
               ESCF=(EE+ENUCLR)*23.061 +ATHEAT
               WRITE(6,'(27X,''AFTER MECI, ENERGY  '',F14.7)')ESCF
            ENDIF
         ENDIF
      ENDIF
      NSCF=NSCF+1
cmgc      TITER2=SECOND()
      IF(TIMES) WRITE(6,'('' TIME FOR SCF CALCULATION'',F8.2,
     1''    INTEGRAL'',F8.2)')TITER2-TITER1,TITER2-TIME0
      IF(DEBUG)WRITE(6,'('' NO. OF ITERATIONS ='',I3)')NITER
C            IF(FORCE)  SCFCRT=1.D-5
      IF(ALLCON.AND.ABS(BSHIFT-4.44d0).LT.1.D-7)THEN
         CAMKIN=.FALSE.
         ALLCON=.FALSE.
         NEWDG=.FALSE.
         BSHIFT=-10.D0
         OKPULY=.FALSE.
      ENDIF
      SHIFT=1.D0
cmgc
         e_cory = escf
cmgc
      RETURN
C
      END
      SUBROUTINE JAB(IA,JA,LLPERM,JINDEX, JJNDEX,PJA,PJB,W, F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION LLPERM(10), PJA(16), PJB(16), W(*), F(*),
     +JINDEX(256), JJNDEX(256)
C    +, SUMA(10), SUMB(10)
         SAVE                                                           GL0892
C
C  FOR VECTOR MACHINES: REMOVE THE ARRAYS SUMA AND SUMB IN DIMENSION STM
C                       UNCOMMENT LINES BETWEEN VECTOR CODE BOUNDARIES
C                       COMMENT LINES BETWEEN SCALAR CODE BOUNDARIES
C  FOR SCALAR MACHINES: ADD THE ARRAYS SUMA AND SUMB IN DIMENSION STMT
C                       UNCOMMENT LINES BETWEEN SCALAR CODE BOUNDARIES
C                       COMMENT LINES BETWEEN VECTOR CODE BOUNDARIES
C
C the following line is the start of the vector code
                  I=0
                  DO 100 I5=1,4
                  IIA=IA+I5-1
                  IJA=JA+I5-1
                  IOFF=(IIA*(IIA-1))/2+IA-1
                  JOFF=(IJA*(IJA-1))/2+JA-1
                  DO 100 I6=1,I5
                  IOFF=IOFF+1
                  JOFF=JOFF+1
                        I=I+1
                        L=LLPERM(I)
                        SUMA=0
                        SUMB=0
                        DO 90 K=1,16
                           L=L+1
                           SUMB=SUMB+PJA(K)*W(JJNDEX(L))
   90                   SUMA=SUMA+PJB(K)*W(JINDEX(L))
                        F(IOFF)=F(IOFF)+SUMA
  100             F(JOFF)=F(JOFF)+SUMB
C the previous line is the end of the vector code
C the following line is the start of the scalar code
C     SUMA( 1)=
C    1+PJA( 1)*W(  1)+PJA( 2)*W( 11)+PJA( 3)*W( 31)+PJA( 4)*W( 61)
C    1+PJA( 5)*W( 11)+PJA( 6)*W( 21)+PJA( 7)*W( 41)+PJA( 8)*W( 71)
C    1+PJA( 9)*W( 31)+PJA(10)*W( 41)+PJA(11)*W( 51)+PJA(12)*W( 81)
C    1+PJA(13)*W( 61)+PJA(14)*W( 71)+PJA(15)*W( 81)+PJA(16)*W( 91)
C     SUMA( 2)=
C    1+PJA( 1)*W(  2)+PJA( 2)*W( 12)+PJA( 3)*W( 32)+PJA( 4)*W( 62)
C    1+PJA( 5)*W( 12)+PJA( 6)*W( 22)+PJA( 7)*W( 42)+PJA( 8)*W( 72)
C    1+PJA( 9)*W( 32)+PJA(10)*W( 42)+PJA(11)*W( 52)+PJA(12)*W( 82)
C    1+PJA(13)*W( 62)+PJA(14)*W( 72)+PJA(15)*W( 82)+PJA(16)*W( 92)
C     SUMA( 3)=
C    1+PJA( 1)*W(  3)+PJA( 2)*W( 13)+PJA( 3)*W( 33)+PJA( 4)*W( 63)
C    1+PJA( 5)*W( 13)+PJA( 6)*W( 23)+PJA( 7)*W( 43)+PJA( 8)*W( 73)
C    1+PJA( 9)*W( 33)+PJA(10)*W( 43)+PJA(11)*W( 53)+PJA(12)*W( 83)
C    1+PJA(13)*W( 63)+PJA(14)*W( 73)+PJA(15)*W( 83)+PJA(16)*W( 93)
C     SUMA( 4)=
C    1+PJA( 1)*W(  4)+PJA( 2)*W( 14)+PJA( 3)*W( 34)+PJA( 4)*W( 64)
C    1+PJA( 5)*W( 14)+PJA( 6)*W( 24)+PJA( 7)*W( 44)+PJA( 8)*W( 74)
C    1+PJA( 9)*W( 34)+PJA(10)*W( 44)+PJA(11)*W( 54)+PJA(12)*W( 84)
C    1+PJA(13)*W( 64)+PJA(14)*W( 74)+PJA(15)*W( 84)+PJA(16)*W( 94)
C     SUMA( 5)=
C    1+PJA( 1)*W(  5)+PJA( 2)*W( 15)+PJA( 3)*W( 35)+PJA( 4)*W( 65)
C    1+PJA( 5)*W( 15)+PJA( 6)*W( 25)+PJA( 7)*W( 45)+PJA( 8)*W( 75)
C    1+PJA( 9)*W( 35)+PJA(10)*W( 45)+PJA(11)*W( 55)+PJA(12)*W( 85)
C    1+PJA(13)*W( 65)+PJA(14)*W( 75)+PJA(15)*W( 85)+PJA(16)*W( 95)
C     SUMA( 6)=
C    1+PJA( 1)*W(  6)+PJA( 2)*W( 16)+PJA( 3)*W( 36)+PJA( 4)*W( 66)
C    1+PJA( 5)*W( 16)+PJA( 6)*W( 26)+PJA( 7)*W( 46)+PJA( 8)*W( 76)
C    1+PJA( 9)*W( 36)+PJA(10)*W( 46)+PJA(11)*W( 56)+PJA(12)*W( 86)
C    1+PJA(13)*W( 66)+PJA(14)*W( 76)+PJA(15)*W( 86)+PJA(16)*W( 96)
C     SUMA( 7)=
C    1+PJA( 1)*W(  7)+PJA( 2)*W( 17)+PJA( 3)*W( 37)+PJA( 4)*W( 67)
C    1+PJA( 5)*W( 17)+PJA( 6)*W( 27)+PJA( 7)*W( 47)+PJA( 8)*W( 77)
C    1+PJA( 9)*W( 37)+PJA(10)*W( 47)+PJA(11)*W( 57)+PJA(12)*W( 87)
C    1+PJA(13)*W( 67)+PJA(14)*W( 77)+PJA(15)*W( 87)+PJA(16)*W( 97)
C     SUMA( 8)=
C    1+PJA( 1)*W(  8)+PJA( 2)*W( 18)+PJA( 3)*W( 38)+PJA( 4)*W( 68)
C    1+PJA( 5)*W( 18)+PJA( 6)*W( 28)+PJA( 7)*W( 48)+PJA( 8)*W( 78)
C    1+PJA( 9)*W( 38)+PJA(10)*W( 48)+PJA(11)*W( 58)+PJA(12)*W( 88)
C    1+PJA(13)*W( 68)+PJA(14)*W( 78)+PJA(15)*W( 88)+PJA(16)*W( 98)
C     SUMA( 9)=
C    1+PJA( 1)*W(  9)+PJA( 2)*W( 19)+PJA( 3)*W( 39)+PJA( 4)*W( 69)
C    1+PJA( 5)*W( 19)+PJA( 6)*W( 29)+PJA( 7)*W( 49)+PJA( 8)*W( 79)
C    1+PJA( 9)*W( 39)+PJA(10)*W( 49)+PJA(11)*W( 59)+PJA(12)*W( 89)
C    1+PJA(13)*W( 69)+PJA(14)*W( 79)+PJA(15)*W( 89)+PJA(16)*W( 99)
C     SUMA(10)=
C    1+PJA( 1)*W( 10)+PJA( 2)*W( 20)+PJA( 3)*W( 40)+PJA( 4)*W( 70)
C    1+PJA( 5)*W( 20)+PJA( 6)*W( 30)+PJA( 7)*W( 50)+PJA( 8)*W( 80)
C    1+PJA( 9)*W( 40)+PJA(10)*W( 50)+PJA(11)*W( 60)+PJA(12)*W( 90)
C    1+PJA(13)*W( 70)+PJA(14)*W( 80)+PJA(15)*W( 90)+PJA(16)*W(100)
C     SUMB( 1)=
C    1+PJB( 1)*W(  1)+PJB( 2)*W(  2)+PJB( 3)*W(  4)+PJB( 4)*W(  7)
C    1+PJB( 5)*W(  2)+PJB( 6)*W(  3)+PJB( 7)*W(  5)+PJB( 8)*W(  8)
C    1+PJB( 9)*W(  4)+PJB(10)*W(  5)+PJB(11)*W(  6)+PJB(12)*W(  9)
C    1+PJB(13)*W(  7)+PJB(14)*W(  8)+PJB(15)*W(  9)+PJB(16)*W( 10)
C     SUMB( 2)=
C    1+PJB( 1)*W( 11)+PJB( 2)*W( 12)+PJB( 3)*W( 14)+PJB( 4)*W( 17)
C    1+PJB( 5)*W( 12)+PJB( 6)*W( 13)+PJB( 7)*W( 15)+PJB( 8)*W( 18)
C    1+PJB( 9)*W( 14)+PJB(10)*W( 15)+PJB(11)*W( 16)+PJB(12)*W( 19)
C    1+PJB(13)*W( 17)+PJB(14)*W( 18)+PJB(15)*W( 19)+PJB(16)*W( 20)
C     SUMB( 3)=
C    1+PJB( 1)*W( 21)+PJB( 2)*W( 22)+PJB( 3)*W( 24)+PJB( 4)*W( 27)
C    1+PJB( 5)*W( 22)+PJB( 6)*W( 23)+PJB( 7)*W( 25)+PJB( 8)*W( 28)
C    1+PJB( 9)*W( 24)+PJB(10)*W( 25)+PJB(11)*W( 26)+PJB(12)*W( 29)
C    1+PJB(13)*W( 27)+PJB(14)*W( 28)+PJB(15)*W( 29)+PJB(16)*W( 30)
C     SUMB( 4)=
C    1+PJB( 1)*W( 31)+PJB( 2)*W( 32)+PJB( 3)*W( 34)+PJB( 4)*W( 37)
C    1+PJB( 5)*W( 32)+PJB( 6)*W( 33)+PJB( 7)*W( 35)+PJB( 8)*W( 38)
C    1+PJB( 9)*W( 34)+PJB(10)*W( 35)+PJB(11)*W( 36)+PJB(12)*W( 39)
C    1+PJB(13)*W( 37)+PJB(14)*W( 38)+PJB(15)*W( 39)+PJB(16)*W( 40)
C     SUMB( 5)=
C    1+PJB( 1)*W( 41)+PJB( 2)*W( 42)+PJB( 3)*W( 44)+PJB( 4)*W( 47)
C    1+PJB( 5)*W( 42)+PJB( 6)*W( 43)+PJB( 7)*W( 45)+PJB( 8)*W( 48)
C    1+PJB( 9)*W( 44)+PJB(10)*W( 45)+PJB(11)*W( 46)+PJB(12)*W( 49)
C    1+PJB(13)*W( 47)+PJB(14)*W( 48)+PJB(15)*W( 49)+PJB(16)*W( 50)
C     SUMB( 6)=
C    1+PJB( 1)*W( 51)+PJB( 2)*W( 52)+PJB( 3)*W( 54)+PJB( 4)*W( 57)
C    1+PJB( 5)*W( 52)+PJB( 6)*W( 53)+PJB( 7)*W( 55)+PJB( 8)*W( 58)
C    1+PJB( 9)*W( 54)+PJB(10)*W( 55)+PJB(11)*W( 56)+PJB(12)*W( 59)
C    1+PJB(13)*W( 57)+PJB(14)*W( 58)+PJB(15)*W( 59)+PJB(16)*W( 60)
C     SUMB( 7)=
C    1+PJB( 1)*W( 61)+PJB( 2)*W( 62)+PJB( 3)*W( 64)+PJB( 4)*W( 67)
C    1+PJB( 5)*W( 62)+PJB( 6)*W( 63)+PJB( 7)*W( 65)+PJB( 8)*W( 68)
C    1+PJB( 9)*W( 64)+PJB(10)*W( 65)+PJB(11)*W( 66)+PJB(12)*W( 69)
C    1+PJB(13)*W( 67)+PJB(14)*W( 68)+PJB(15)*W( 69)+PJB(16)*W( 70)
C     SUMB( 8)=
C    1+PJB( 1)*W( 71)+PJB( 2)*W( 72)+PJB( 3)*W( 74)+PJB( 4)*W( 77)
C    1+PJB( 5)*W( 72)+PJB( 6)*W( 73)+PJB( 7)*W( 75)+PJB( 8)*W( 78)
C    1+PJB( 9)*W( 74)+PJB(10)*W( 75)+PJB(11)*W( 76)+PJB(12)*W( 79)
C    1+PJB(13)*W( 77)+PJB(14)*W( 78)+PJB(15)*W( 79)+PJB(16)*W( 80)
C     SUMB( 9)=
C    1+PJB( 1)*W( 81)+PJB( 2)*W( 82)+PJB( 3)*W( 84)+PJB( 4)*W( 87)
C    1+PJB( 5)*W( 82)+PJB( 6)*W( 83)+PJB( 7)*W( 85)+PJB( 8)*W( 88)
C    1+PJB( 9)*W( 84)+PJB(10)*W( 85)+PJB(11)*W( 86)+PJB(12)*W( 89)
C    1+PJB(13)*W( 87)+PJB(14)*W( 88)+PJB(15)*W( 89)+PJB(16)*W( 90)
C     SUMB(10)=
C    1+PJB( 1)*W( 91)+PJB( 2)*W( 92)+PJB( 3)*W( 94)+PJB( 4)*W( 97)
C    1+PJB( 5)*W( 92)+PJB( 6)*W( 93)+PJB( 7)*W( 95)+PJB( 8)*W( 98)
C    1+PJB( 9)*W( 94)+PJB(10)*W( 95)+PJB(11)*W( 96)+PJB(12)*W( 99)
C    1+PJB(13)*W( 97)+PJB(14)*W( 98)+PJB(15)*W( 99)+PJB(16)*W(100)
C                 I=0
C                 DO 100 I5=1,4
C                 IIA=IA+I5-1
C                 IJA=JA+I5-1
C                 IOFF=(IIA*(IIA-1))/2+IA-1
C                 JOFF=(IJA*(IJA-1))/2+JA-1
C                 DO 100 I6=1,I5
C                 IOFF=IOFF+1
C                 JOFF=JOFF+1
C                       I=I+1
C                       F(IOFF)=F(IOFF)+SUMB(I)
C 100             F(JOFF)=F(JOFF)+SUMA(I)
C the previous line is the end of the scalar code
      RETURN
      END
      SUBROUTINE KAB(IA,JA, PK, W, KINDEX, F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PK(*), W(*), F(*),KINDEX(256)
C    +, SUM(16)
         SAVE                                                           GL0892
C
C  FOR VECTOR MACHINES: REMOVE THE ARRAY SUM IN DIMENSION STMT
C                       UNCOMMENT LINES BETWEEN VECTOR CODE BOUNDARIES
C                       COMMENT LINES BETWEEN SCALAR CODE BOUNDARIES
C  FOR SCALAR MACHINES: ADD THE ARRAY SUM IN DIMENSION STMT
C                       UNCOMMENT LINES BETWEEN SCALAR CODE BOUNDARIES
C                       COMMENT LINES BETWEEN VECTOR CODE BOUNDARIES
C
C the following line is the start of the vector code
                  L=0
                  M=0
                  DO 130 J1=IA,IA+3
                  J=(J1*(J1-1))/2
                  DO 130 J2=JA,JA+3
                  M=M+1
                  J3=J+J2
                     SUM=0
                     DO 120 I=1,16
                        L=L+1
  120                SUM=SUM+PK(I)*W(KINDEX(L))
  130             F(J3)=F(J3)-SUM
C the previous line is the end of the vector code
C the following line is the start of the scalar code
C     SUM( 1)=
C    1+PK( 1)*W(  1)+PK( 2)*W(  2)+PK( 3)*W(  4)+PK( 4)*W(  7)
C    1+PK( 5)*W( 11)+PK( 6)*W( 12)+PK( 7)*W( 14)+PK( 8)*W( 17)
C    1+PK( 9)*W( 31)+PK(10)*W( 32)+PK(11)*W( 34)+PK(12)*W( 37)
C    1+PK(13)*W( 61)+PK(14)*W( 62)+PK(15)*W( 64)+PK(16)*W( 67)
C     SUM( 2)=
C    1+PK( 1)*W(  2)+PK( 2)*W(  3)+PK( 3)*W(  5)+PK( 4)*W(  8)
C    1+PK( 5)*W( 12)+PK( 6)*W( 13)+PK( 7)*W( 15)+PK( 8)*W( 18)
C    1+PK( 9)*W( 32)+PK(10)*W( 33)+PK(11)*W( 35)+PK(12)*W( 38)
C    1+PK(13)*W( 62)+PK(14)*W( 63)+PK(15)*W( 65)+PK(16)*W( 68)
C     SUM( 3)=
C    1+PK( 1)*W(  4)+PK( 2)*W(  5)+PK( 3)*W(  6)+PK( 4)*W(  9)
C    1+PK( 5)*W( 14)+PK( 6)*W( 15)+PK( 7)*W( 16)+PK( 8)*W( 19)
C    1+PK( 9)*W( 34)+PK(10)*W( 35)+PK(11)*W( 36)+PK(12)*W( 39)
C    1+PK(13)*W( 64)+PK(14)*W( 65)+PK(15)*W( 66)+PK(16)*W( 69)
C     SUM( 4)=
C    1+PK( 1)*W(  7)+PK( 2)*W(  8)+PK( 3)*W(  9)+PK( 4)*W( 10)
C    1+PK( 5)*W( 17)+PK( 6)*W( 18)+PK( 7)*W( 19)+PK( 8)*W( 20)
C    1+PK( 9)*W( 37)+PK(10)*W( 38)+PK(11)*W( 39)+PK(12)*W( 40)
C    1+PK(13)*W( 67)+PK(14)*W( 68)+PK(15)*W( 69)+PK(16)*W( 70)
C     SUM( 5)=
C    1+PK( 1)*W( 11)+PK( 2)*W( 12)+PK( 3)*W( 14)+PK( 4)*W( 17)
C    1+PK( 5)*W( 21)+PK( 6)*W( 22)+PK( 7)*W( 24)+PK( 8)*W( 27)
C    1+PK( 9)*W( 41)+PK(10)*W( 42)+PK(11)*W( 44)+PK(12)*W( 47)
C    1+PK(13)*W( 71)+PK(14)*W( 72)+PK(15)*W( 74)+PK(16)*W( 77)
C     SUM( 6)=
C    1+PK( 1)*W( 12)+PK( 2)*W( 13)+PK( 3)*W( 15)+PK( 4)*W( 18)
C    1+PK( 5)*W( 22)+PK( 6)*W( 23)+PK( 7)*W( 25)+PK( 8)*W( 28)
C    1+PK( 9)*W( 42)+PK(10)*W( 43)+PK(11)*W( 45)+PK(12)*W( 48)
C    1+PK(13)*W( 72)+PK(14)*W( 73)+PK(15)*W( 75)+PK(16)*W( 78)
C     SUM( 7)=
C    1+PK( 1)*W( 14)+PK( 2)*W( 15)+PK( 3)*W( 16)+PK( 4)*W( 19)
C    1+PK( 5)*W( 24)+PK( 6)*W( 25)+PK( 7)*W( 26)+PK( 8)*W( 29)
C    1+PK( 9)*W( 44)+PK(10)*W( 45)+PK(11)*W( 46)+PK(12)*W( 49)
C    1+PK(13)*W( 74)+PK(14)*W( 75)+PK(15)*W( 76)+PK(16)*W( 79)
C     SUM( 8)=
C    1+PK( 1)*W( 17)+PK( 2)*W( 18)+PK( 3)*W( 19)+PK( 4)*W( 20)
C    1+PK( 5)*W( 27)+PK( 6)*W( 28)+PK( 7)*W( 29)+PK( 8)*W( 30)
C    1+PK( 9)*W( 47)+PK(10)*W( 48)+PK(11)*W( 49)+PK(12)*W( 50)
C    1+PK(13)*W( 77)+PK(14)*W( 78)+PK(15)*W( 79)+PK(16)*W( 80)
C     SUM( 9)=
C    1+PK( 1)*W( 31)+PK( 2)*W( 32)+PK( 3)*W( 34)+PK( 4)*W( 37)
C    1+PK( 5)*W( 41)+PK( 6)*W( 42)+PK( 7)*W( 44)+PK( 8)*W( 47)
C    1+PK( 9)*W( 51)+PK(10)*W( 52)+PK(11)*W( 54)+PK(12)*W( 57)
C    1+PK(13)*W( 81)+PK(14)*W( 82)+PK(15)*W( 84)+PK(16)*W( 87)
C     SUM(10)=
C    1+PK( 1)*W( 32)+PK( 2)*W( 33)+PK( 3)*W( 35)+PK( 4)*W( 38)
C    1+PK( 5)*W( 42)+PK( 6)*W( 43)+PK( 7)*W( 45)+PK( 8)*W( 48)
C    1+PK( 9)*W( 52)+PK(10)*W( 53)+PK(11)*W( 55)+PK(12)*W( 58)
C    1+PK(13)*W( 82)+PK(14)*W( 83)+PK(15)*W( 85)+PK(16)*W( 88)
C     SUM(11)=
C    1+PK( 1)*W( 34)+PK( 2)*W( 35)+PK( 3)*W( 36)+PK( 4)*W( 39)
C    1+PK( 5)*W( 44)+PK( 6)*W( 45)+PK( 7)*W( 46)+PK( 8)*W( 49)
C    1+PK( 9)*W( 54)+PK(10)*W( 55)+PK(11)*W( 56)+PK(12)*W( 59)
C    1+PK(13)*W( 84)+PK(14)*W( 85)+PK(15)*W( 86)+PK(16)*W( 89)
C     SUM(12)=
C    1+PK( 1)*W( 37)+PK( 2)*W( 38)+PK( 3)*W( 39)+PK( 4)*W( 40)
C    1+PK( 5)*W( 47)+PK( 6)*W( 48)+PK( 7)*W( 49)+PK( 8)*W( 50)
C    1+PK( 9)*W( 57)+PK(10)*W( 58)+PK(11)*W( 59)+PK(12)*W( 60)
C    1+PK(13)*W( 87)+PK(14)*W( 88)+PK(15)*W( 89)+PK(16)*W( 90)
C     SUM(13)=
C    1+PK( 1)*W( 61)+PK( 2)*W( 62)+PK( 3)*W( 64)+PK( 4)*W( 67)
C    1+PK( 5)*W( 71)+PK( 6)*W( 72)+PK( 7)*W( 74)+PK( 8)*W( 77)
C    1+PK( 9)*W( 81)+PK(10)*W( 82)+PK(11)*W( 84)+PK(12)*W( 87)
C    1+PK(13)*W( 91)+PK(14)*W( 92)+PK(15)*W( 94)+PK(16)*W( 97)
C     SUM(14)=
C    1+PK( 1)*W( 62)+PK( 2)*W( 63)+PK( 3)*W( 65)+PK( 4)*W( 68)
C    1+PK( 5)*W( 72)+PK( 6)*W( 73)+PK( 7)*W( 75)+PK( 8)*W( 78)
C    1+PK( 9)*W( 82)+PK(10)*W( 83)+PK(11)*W( 85)+PK(12)*W( 88)
C    1+PK(13)*W( 92)+PK(14)*W( 93)+PK(15)*W( 95)+PK(16)*W( 98)
C     SUM(15)=
C    1+PK( 1)*W( 64)+PK( 2)*W( 65)+PK( 3)*W( 66)+PK( 4)*W( 69)
C    1+PK( 5)*W( 74)+PK( 6)*W( 75)+PK( 7)*W( 76)+PK( 8)*W( 79)
C    1+PK( 9)*W( 84)+PK(10)*W( 85)+PK(11)*W( 86)+PK(12)*W( 89)
C    1+PK(13)*W( 94)+PK(14)*W( 95)+PK(15)*W( 96)+PK(16)*W( 99)
C     SUM(16)=
C    1+PK( 1)*W( 67)+PK( 2)*W( 68)+PK( 3)*W( 69)+PK( 4)*W( 70)
C    1+PK( 5)*W( 77)+PK( 6)*W( 78)+PK( 7)*W( 79)+PK( 8)*W( 80)
C    1+PK( 9)*W( 87)+PK(10)*W( 88)+PK(11)*W( 89)+PK(12)*W( 90)
C    1+PK(13)*W( 97)+PK(14)*W( 98)+PK(15)*W( 99)+PK(16)*W(100)
C                  M=0
C                 DO 130 J1=IA,IA+3
C                 J=(J1*(J1-1))/2
C                 DO 130 J2=JA,JA+3
C                 M=M+1
C                 J3=J+J2
C 130             F(J3)=F(J3)-SUM(M)
C the previous line is the end of the scalar code
       RETURN
      END
      SUBROUTINE LINMIN(XPARAM,STEP,PVECT,NVAR,FUNCT,OKF,OKC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION XPARAM(NVAR),PVECT(NVAR)
      COMMON /GRAVEC/ COSINE
      COMMON /NUMCAL/ NUMCAL
C*********************************************************************
C
C  LINMIN DOES A LINE MINIMISATION.
C
C  ON INPUT:  XPARAM = STARTING COORDINATE OF SEARCH.
C             STEP   = STEP SIZE FOR INITIATING SEARCH.
C             PVECT  = DIRECTION OF SEARCH.
C             NVAR   = NUMBER OF VARIABLES IN XPARAM.
C             FUNCT  = INITIAL VALUE OF THE FUNCTION TO BE MINIMIZED.
C             ISOK   = NOT IMPORTANT.
C             COSINE = COSINE OF ANGLE OF CURRENT AND PREVIOUS GRADIENT.
C
C  ON OUTPUT: XPARAM = COORDINATE OF MINIMUM OF FUNCTI0N.
C             STEP   = NEW STEP SIZE, USED IN NEXT CALL OF LINMIN.
C             PVECT  = UNCHANGED, OR NEGATED, DEPENDING ON STEP.
C             FUNCT  = FINAL, MINIMUM VALUE OF THE FUNCTION.
C             OKF    = TRUE IF LINMIN IMPROVED FUNCT, FALSE OTHERWISE.
C             OKC    = TRUE IF LINMIN FOUND THE MINIMUM, FALSE OTHERWISE
C
C**********************************************************************
      COMMON /KEYWRD/ KEYWRD
C
C  THE FOLLOWING COMMON IS USED TO FIND OUT IF A NON-VARIATIONALLY
C  OPTIMIZED WAVE-FUNCTION IS BEING USED.
C
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      CHARACTER*80 KEYWRD
      DIMENSION PHI(3), VT(3)
      DIMENSION XSTOR(MAXPAR)
      INTEGER LEFT,RIGHT,CENTER
      LOGICAL PRINT,OKF,OKC, FULSCF, ASKFUL, HALFE
         SAVE                                                           GL0892
      DATA ICALCN /0/
      IF (ICALCN.NE.NUMCAL) THEN
         HALFE=(INDEX(KEYWRD,'C.I.') .NE. 0 .OR. NCLOSE.NE.NOPEN)
         ASKFUL=(INDEX(KEYWRD,'FULSCF') .NE. 0)
         DROP=0.00002D0
         IF(INDEX(KEYWRD,'PREC') .NE. 0) DROP=DROP*0.01D0
         IF(INDEX(KEYWRD,'GNORM') .NE. 0)
     1DROP=DROP*MIN(READA(KEYWRD,INDEX(KEYWRD,'GNORM')),1.D0)
         XMAXM  = 0.4D0
         I      = 2
         STEP   = 1.D0
         MAXLIN = 15
         XCRIT  = 0.0001D0
         IF(INDEX(KEYWRD,'FORCE') .NE. 0) THEN
            I=3
            XCRIT=0.00001D0
         ENDIF
         ANGLE=0.8D0
         IF(HALFE) ANGLE=-2.D0
         COSINE=99.99D0
C
C  ANGLE IS USED TO DECIDE IF P IS TO BE UPDATED AS CALCULATION
C        PROCEEDS.
C
         YMAXST  = 0.4D0
         EPS=10**(-I)
         TEE=EPS
         PRINT=(INDEX(KEYWRD,'LINMIN') .NE. 0)
         ICALCN=NUMCAL
      ENDIF
      FULSCF=(ASKFUL.OR.COSINE.GT.ANGLE)
      IF(PRINT) WRITE(6,'(''  FULL SCF CALCULATIONS:'',3X,L1)')FULSCF
      XMAXM=0.D0
      DO 10 I=1,NVAR
         PABS=ABS(PVECT(I))
   10 XMAXM=MAX(XMAXM,PABS)
      XMINM=XMAXM
      XMAXM=YMAXST/XMAXM
      FIN=FUNCT
      SSQLST=FUNCT
      IQUIT=0
      PHI(1)=FUNCT
      VT(1)=0.0D00
      VT(2)=STEP/4.0D00
      IF (VT(2).GT.XMAXM) VT(2)=XMAXM
      FMAX=FUNCT
      FMIN=FUNCT
      STEP=VT(2)
      DO 20 I=1,NVAR
   20 XPARAM(I)=XPARAM(I)+STEP*PVECT(I)
      CALL COMPFG(XPARAM, .TRUE., PHI(2),FULSCF,GRAD,.FALSE.)
      IF(PHI(2).GT.FMAX) FMAX=PHI(2)
      IF(PHI(2).LT.FMIN) FMIN=PHI(2)
      sqstor = phi(2)
      estor = energy
      alfs = step
      do 25 iexch = 1,nvar
      xstor(iexch) = xparam(iexch)
25    continue
      IF (PHI(1).LE.PHI(2)) GO TO 30
      GO TO 40
   30 VT(3)=-VT(2)
      LEFT=3
      CENTER=1
      RIGHT=2
      GO TO 50
   40 VT(3)=2.0D00*VT(2)
      LEFT=1
      CENTER=2
      RIGHT=3
   50 STLAST=VT(3)
      STEP=STLAST-STEP
      DO 60 I=1,NVAR
   60 XPARAM(I)=XPARAM(I)+STEP*PVECT(I)
      CALL COMPFG (XPARAM, .TRUE., FUNCT,FULSCF,GRAD,.FALSE.)
      IF(FUNCT.GT.FMAX) FMAX=FUNCT
      IF(FUNCT.LT.FMIN) FMIN=FUNCT
      if(funct .lt. sqstor) then
        sqstor = funct
        estor = energy
        alfs = step
        do 65 iexch = 1,nvar
        xstor(iexch) = xparam(iexch)
65      continue
      endif
      IF (FUNCT.LT.FIN) IQUIT=1
      PHI(3)=FUNCT
      IF (PRINT)WRITE (6,230) VT(1),PHI(1),VT(2),PHI(2),VT(3),PHI(3)
      OKC=.TRUE.
      DO 180 ICTR=3,MAXLIN
         ALPHA=VT(2)-VT(3)
         BETA=VT(3)-VT(1)
         GAMMA=VT(1)-VT(2)
         IF(ABS(ALPHA*BETA*GAMMA) .GT. 1.D-20)THEN
            ALPHA=-(PHI(1)*ALPHA+PHI(2)*BETA+PHI(3)*GAMMA)/(ALPHA*BETA*G
     1AMM   A)
         ELSE
            GOTO 190
         ENDIF
         BETA=((PHI(1)-PHI(2))/GAMMA)-ALPHA*(VT(1)+VT(2))
         IF (ALPHA) 70,70,100
   70    IF (PHI(RIGHT).GT.PHI(LEFT)) GO TO 80
         STEP=3.0D00*VT(RIGHT)-2.0D00*VT(CENTER)
         GO TO 90
   80    STEP=3.0D00*VT(LEFT)-2.0D00*VT(CENTER)
   90    S=STEP-STLAST
         IF (ABS(S).GT.XMAXM) S=SIGN(XMAXM,S)*(1+0.01*(XMAXM/S))
         STEP=S+STLAST
         GO TO 110
  100    STEP=-BETA/(2.0D00*ALPHA)
         S=STEP-STLAST
         XXM=2.0D00*XMAXM
         IF (ABS(S).GT.XXM) S=SIGN(XXM,S)*(1+0.01*(XXM/S))
         STEP=S+STLAST
  110    CONTINUE
         IF (ICTR.LE.3) GO TO 120
         AABS=ABS(S*XMINM)
         IF (AABS.LT.XCRIT) GO TO 190
  120    CONTINUE
         DO 130 I=1,NVAR
  130    XPARAM(I)=XPARAM(I)+S*PVECT(I)
         FUNOLD=FUNCT
         CALL COMPFG (XPARAM, .TRUE., FUNCT,FULSCF,GRAD,.FALSE.)
         IF(FUNCT.GT.FMAX) FMAX=FUNCT
         IF(FUNCT.LT.FMIN) FMIN=FUNCT
      if(funct .lt. sqstor) then
        sqstor = funct
        estor = energy
        alfs = step
        do 135 iexch = 1,nvar
        xstor(iexch) = xparam(iexch)
135     continue
      endif
         IF (FUNCT.LT.FIN) IQUIT=1
         IF (PRINT) WRITE (6,240) VT(LEFT),PHI(LEFT),
     1                            VT(CENTER),PHI(CENTER),
     2                            VT(RIGHT),PHI(RIGHT),STEP,FUNCT
C
C TEST TO EXIT FROM LINMIN IF NOT DROPPING IN VALUE OF FUNCTION FAST.
C
         TINY = MAX((SSQLST-FMIN)*0.2D0 , DROP)
         TINY = MIN(TINY,0.5D0)
         IF(PRINT) WRITE(6,'(''  TINY'',F14.9)')TINY
         IF(ABS(FUNOLD-FUNCT) .LT. TINY .AND. IQUIT .EQ. 1) GOTO 190
         IF ((ABS(STEP-STLAST).LE.EPS*ABS(STEP+STLAST)+TEE).
     1   AND.(IQUIT.EQ.1)) GO TO 190
         STLAST=STEP
         IF ((STEP.GT.VT(RIGHT)).OR.(STEP.GT.VT(CENTER)
     1        .AND.FUNCT.LT.PHI(CENTER)).OR.(STEP.GT.VT(LEFT)
     2        .AND.STEP.LT.VT(CENTER).AND.FUNCT.GT.PHI(CENTER)))
     3         GOTO 140
         VT(RIGHT)=STEP
         PHI(RIGHT)=FUNCT
         GO TO 150
  140    VT(LEFT)=STEP
         PHI(LEFT)=FUNCT
  150    IF (VT(CENTER).LT.VT(RIGHT)) GO TO 160
         I=CENTER
         CENTER=RIGHT
         RIGHT=I
  160    IF (VT(LEFT).LT.VT(CENTER)) GO TO 170
         I=LEFT
         LEFT=CENTER
         CENTER=I
  170    IF (VT(CENTER).LT.VT(RIGHT)) GO TO 180
         I=CENTER
         CENTER=RIGHT
         RIGHT=I
  180 CONTINUE
      OKC=.FALSE.
  190 CONTINUE
      funct = sqstor
      energy = estor
      step = alfs
      do 195 iexch = 1,nvar
      xparam(iexch) = xstor(iexch)
195   continue
      OKF = (FUNCT.LT.SSQLST)
      IF (FUNCT.GE.SSQLST) RETURN
      IF (STEP) 200,220,220
  200 STEP=-STEP
      DO 210 I=1,NVAR
  210 PVECT(I)=-PVECT(I)
  220 CONTINUE
      RETURN
C
  230 FORMAT ( 11H ---QLINMN ,/5X, 'LEFT   ...',F17.8,F17.11/5X,
     1 'CENTER ...',F17.8,F17.11,/5X, 'RIGHT  ...',F17.8,F17.11,/)
  240 FORMAT (5X,'LEFT    ...',F17.8,F17.11,/5X,'CENTER  ...',
     1F17.8,F17.11,/5X,'RIGHT   ...',F17.8,F17.11,/5X,
     2 'NEW     ...',F17.8,F17.11,/)
C
      END
      SUBROUTINE LOCAL(C,MDIM,NOCC,EIG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION C(MDIM,MDIM), EIG(MAXORB)
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
C**********************************************************************
C
C   LOCALISATION SUBROUTINE
C ON INPUT
C        C = EIGENVECTORS IN AN MDIM*MDIM MATRIX
C        NOCC = NUMBER OF FILLED LEVELS
C        NORBS = NUMBER OF ORBITALS
C        NUMAT = NUMBER OF ATOMS
C        NLAST   = INTEGER ARRAY OF ATOM ORBITAL COUNTERS
C        NFIRST   = INTEGER ARRAY OF ATOM ORBITAL COUNTERS
C
C       SUBROUTINE MAXIMIZES (PSI)**4
C       REFERENCE_
C       A NEW RAPID METHOD FOR ORBITAL LOCALISATION, P.G. PERKINS AND
C       J.J.P. STEWART, J.C.S. FARADAY (II) 77, 000, (1981).
C
C       MODIFIED AND CORRECTED TO AVOID SIGMA-PI ORBITAL MIXING BY
C       JUAN CARLOS PANIAGUA, UNIVERSITY OF BARCELONA, MAY 1983.
C
C**********************************************************************
      COMMON /SCRACH/ COLD(MAXORB,MAXORB),XDUMY(MAXPAR**2-MAXORB*MAXORB)
      DIMENSION EIG1(MAXORB),PSI1(MAXORB),PSI2(MAXORB),
     1          CII(MAXORB), REFEIG(MAXORB),IEL(20)
      CHARACTER*2 ELEMNT(99)
         SAVE                                                           GL0892
      DATA ELEMNT/'H','HE',
     1 'LI','BE','B','C','N','O','F','NE',
     2 'NA','MG','AL','SI','P','S','CL','AR',
     3 'K','CA','SC','TI','V','CR','MN','FE','CO','NI','CU',
     4 'ZN','GA','GE','AS','SE','BR','KR',
     5 'RB','SR','Y','ZR','NB','MO','TC','RU','RH','PD','AG',
     6 'CD','IN','SN','SB','TE','I','XE',
     7 'CS','BA','LA','CE','PR','ND','PM','SM','EU','GD','TB','DY',
     8 'HO','ER','TM','YB','LU','HF','TA','W','RE','OS','IR','PT',
     9 'AU','HG','TL','PB','BI','PO','AT','RN',
     1 'FR','RA','AC','TH','PA','U','NP','PU','AM','CM','BK','CF','XX'/
      NITER=100
      EPS=1.0D-7
      DO 10 I=1,NORBS
         REFEIG(I)=EIG(I)
         DO 10 J=1,NORBS
   10 COLD(I,J)=C(I,J)
      ITER=0
   20 CONTINUE
      SUM=0.D0
      ITER=ITER+1
      DO 80 I=1,NOCC
         DO 70 J=1,NOCC
            IF(J.EQ.I) GOTO 70
            XIJJJ=0.0D0
            XJIII=0.0D0
            XIIII=0.0D0
            XJJJJ=0.0D0
            XIJIJ=0.0D0
            XIIJJ=0.0D0
            DO 30 K=1,NORBS
               PSI1(K)=C(K,I)
   30       PSI2(K)=C(K,J)
C NOW FOLLOWS THE RATE-DETERMINING STEP FOR THE CALCULATION
            DO 50 K1=1,NUMAT
               KL=NFIRST(K1)
               KU=NLAST(K1)
               DIJ=0.D0
               DII=0.D0
               DJJ=0.D0
               DO 40 K=KL,KU
                  DIJ=DIJ+PSI1(K)*PSI2(K)
                  DII=DII+PSI1(K)*PSI1(K)
                  DJJ=DJJ+PSI2(K)*PSI2(K)
   40          CONTINUE
               XIJJJ=XIJJJ+DIJ*DJJ
               XJIII=XJIII+DIJ*DII
               XIIII=XIIII+DII*DII
               XJJJJ=XJJJJ+DJJ*DJJ
               XIJIJ=XIJIJ+DIJ*DIJ
               XIIJJ=XIIJJ+DII*DJJ
   50       CONTINUE
            AIJ=XIJIJ-(XIIII+XJJJJ-2.0D0*XIIJJ)/4.0D0
            BIJ=XJIII-XIJJJ
            CA=SQRT(AIJ*AIJ+BIJ*BIJ)
            SA=AIJ+CA
            IF(SA.LT.1.0D-14) GO TO 70
            SUM=SUM+SA
            CA=-AIJ/CA
            CA=(1.0D0+SQRT((1.0D0+CA)/2.0D0))/2.0D0
            IF((2.0D0*CA-1.0D0)*BIJ.LT.0.0D0)CA=1.0D0-CA
            SA=SQRT(1.0D0-CA)
            CA=SQRT(CA)
            DO 60 K=1,NORBS
               C(K,I)=CA*PSI1(K)+SA*PSI2(K)
   60       C(K,J)=-SA*PSI1(K)+CA*PSI2(K)
   70    CONTINUE
   80 CONTINUE
      SUM1=0.D0
      DO 100 I=1,NOCC
         DO 100 J=1,NUMAT
            IL=NFIRST(J)
            IU=NLAST(J)
            X=0.0
            DO 90 K=IL,IU
   90       X=X+C(K,I)**2
  100 SUM1=SUM1+X*X
      IF(SUM.GT.EPS.AND.ITER.LT.NITER) GO TO 20
      WRITE(6,110)ITER,SUM1
  110 FORMAT(/10X,'NUMBER OF ITERATIONS =',I4/
     110X,'LOCALISATION VALUE =',F14.9,/)
      WRITE(6,120)
  120 FORMAT(3X,'NUMBER OF CENTERS',14X,'(COMPOSITION OF ORBITALS)'//)
      DO 150 I=1,NOCC
         SUM=0.D0
         DO 140 J=1,NOCC
            CO=0.D0
            DO 130 K=1,NORBS
  130       CO=CO+COLD(K,J)*C(K,I)
  140    SUM=SUM+CO*CO*EIG(J)
  150 EIG1(I)=SUM
      DO 180 I=1,NOCC
         X=100.D0
         DO 160 J=I,NOCC
            IF (X.LT.EIG1(J))  GOTO  160
            X=EIG1(J)
            I1=J
  160    CONTINUE
         EIG(I)=EIG1(I1)
         X=EIG1(I1)
         EIG1(I1)=EIG1(I)
         EIG1(I)=X
         DO 170 J=1,NORBS
            X=C(J,I1)
            C(J,I1)=C(J,I)
  170    C(J,I)=X
  180 CONTINUE
      DO 250 I=1,NOCC
         X=0.D0
         DO 200 K1=1,NUMAT
            KL=NFIRST(K1)
            KU=NLAST(K1)
            DII=0.D0
            DO 190 K=KL,KU
  190       DII=DII+C(K,I)**2
            X=X+DII*DII
  200    PSI1(K1)=DII*100.D0
         X=1/X
         DO 220 II=1,NUMAT
            SUM=0.D0
            DO 210 J=1,NUMAT
               IF(PSI1(J).LT.SUM) GOTO 210
               SUM=PSI1(J)
               K=J
  210       CONTINUE
            PSI1(K)=0.D0
            CII(II)=SUM
            IEL(II)=K
            IF(SUM.LT.1.D0) GOTO 230
  220    CONTINUE
  230    CONTINUE
         II=II-1
         WRITE(6,240)X,(ELEMNT(NAT(IEL(K))),IEL(K),CII(K),K=1,II)
  240    FORMAT(F10.4,4(5(3X,A2,I3,F6.2),/10X))
  250 CONTINUE
  260 FORMAT(//20X,20H LOCALIZED ORBITALS   ,//)
      WRITE(6,260)
      CALL MATOUT(C,EIG,NOCC,NORBS,MDIM)
  270 FORMAT(10F12.6)
      DO 280 I=1,NOCC
         EIG(I)=REFEIG(I)
         DO 280 J=1,NORBS
  280 C(J,I)=COLD(J,I)
      RETURN
      END
      SUBROUTINE MAMULT(A,B,C,N,ONE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*),B(*),C(*)
         SAVE                                                           GL0892
************************************************************************
*
*   MAMULT MULTIPLIES A BY B AND PUTS THE RESULT IN C
*
************************************************************************
      L=0
      DO 40 I=1,N
         II=((I-1)*I)/2
         DO 40 J=1,I
            JJ=((J-1)*J)/2
            L=L+1
            SUM=0.D0
            DO 10 K=1,J
   10       SUM=SUM+A(II+K)*B(JJ+K)
            DO 20 K=J+1,I
   20       SUM=SUM+A(II+K)*B(((K-1)*K)/2+J)
            DO 30 K=I+1,N
               KK=(K*(K-1))/2
   30       SUM=SUM+A(KK+I)*B(KK+J)
   40 C(L)=SUM+ONE*C(L)
      RETURN
      END
      SUBROUTINE MATOUT (A,B,NC,NR,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION A(NDIM,NDIM), B(NDIM)
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /ELEMTS/ ELEMNT(107)
cmgc
      common /cory/ e_cory
      COMMON /GRADNT/ GRAD(MAXPAR),GNORM
      logical srp_test
cmgc
C**********************************************************************
C
C      MATOUT PRINTS A SQUARE MATRIX OF EIGENVECTORS AND EIGENVALUES
C
C    ON INPUT A CONTAINS THE MATRIX TO BE PRINTED.
C             B CONTAINS THE EIGENVALUES.
C             NC NUMBER OF MOLECULAR ORBITALS TO BE PRINTED.
C             NR IS THE SIZE OF THE SQUARE ARRAY TO BE PRINTED.
C             NDIM IS THE ACTUAL SIZE OF THE SQUARE ARRAY "A".
C             NFIRST AND NLAST CONTAIN ATOM ORBITAL COUNTERS.
C             NAT = ARRAY OF ATOMIC NUMBERS OF ATOMS.
C
C
C***********************************************************************
      CHARACTER*2 ELEMNT, ATORBS(9), ITEXT(MAXORB), JTEXT(MAXORB)
      DIMENSION NATOM(MAXORB)
         SAVE                                                           GL0892
      DATA ATORBS/' S','PX','PY','PZ','X2','XZ','Z2','YZ','XY'/
      IF(NUMAT.EQ.0)GOTO 30
      IF(NLAST(NUMAT).NE.NR) GOTO 30
      DO 20 I=1,NUMAT
         JLO=NFIRST(I)
         JHI=NLAST(I)
         L=NAT(I)
         K=0
         DO 10 J=JLO,JHI
            K=K+1
            ITEXT(J)=ATORBS(K)
            JTEXT(J)=ELEMNT(L)
            NATOM(J)=I
   10    CONTINUE
   20 CONTINUE
      GOTO 50
   30 CONTINUE
      NR=ABS(NR)
      DO 40 I=1,NR
         ITEXT(I)='  '
         JTEXT(I)='  '
   40 NATOM(I)=I
   50 CONTINUE
      KA=1
      KC=6
   60 KB=MIN0(KC,NC)
      WRITE (6,100) (I,I=KA,KB)
      IF(B(1).NE.0.D0)WRITE (6,110) (B(I),I=KA,KB)
cmgc
c      inquire(file='SRP_STUFF', exist=srp_test)
c      if(srp_test) then
c         open(unit=97,form='unformatted',status='old',file='SRP_STUFF')
c         read(97)
c      else
      open(unit=97,form='unformatted',status='unknown',file='SRP_STUFF')
      rewind 97
c      endif

cmgc      write(97)e_cory,gnorm,(b(i), i=1,7)
      write(97)gnorm,(b(i), i=1,7)

      close(unit=97,status='keep')
      call exit(0)
cmgc
      WRITE (6,120)
      LA=1
      LC=40
   70 LB=MIN0(LC,NR)
      DO 80 I=LA,LB
         IF(ITEXT(I).EQ.' S')WRITE(6,120)
         WRITE (6,130) ITEXT(I),JTEXT(I),NATOM(I),(A(I,J),J=KA,KB)
   80 CONTINUE
      IF (LB.EQ.NR) GO TO 90
      LA=LC+1
      LC=LC+40
      WRITE (6,140)
      GO TO 70
   90 IF (KB.EQ.NC) RETURN
      KA=KC+1
      KC=KC+6
      IF (NR.GT.25) WRITE (6,140)
      GO TO 60
C
  100 FORMAT (////,3X,9H ROOT NO.,I5,9I12)
  110 FORMAT (/8X,10F12.5)
  120 FORMAT (2H  )
  130 FORMAT (2(1X,A2),I4,F10.5,10F12.5)
  140 FORMAT (1H1)
C
      END
      SUBROUTINE MOLVAL(C,NDUMMY,P,NOCC,RHFUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION C(NORBS,NORBS), P(*)
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      DIMENSION VAL(MAXORB)
         SAVE                                                           GL0892
      DO 40 I=1,NOCC
         SUM=0.D0
         DO 30 JJ=1,NUMAT
            JL=NFIRST(JJ)
            JU=NLAST(JJ)
            DO 30 J=JL,JU
               DO 30 KK=1,NUMAT
                  IF(KK.EQ.JJ) GOTO 20
                  KL=NFIRST(KK)
                  KU=NLAST(KK)
                  DO 10 K=KL,KU
                     L1=MAX(J,K)
                     L2=J+K-L1
                     L=(L1*(L1-1))/2+L2
                     SUM=SUM+C(J,I)*C(K,I)*P(L)
   10             CONTINUE
   20             CONTINUE
   30    CONTINUE
         VAL(I)=SUM*RHFUHF
   40 CONTINUE
      WRITE(6,'(10F8.4)')(VAL(I),I=1,NOCC)
      END
      SUBROUTINE NUCHAR(LINE,VALUE,NVALUE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
************************************************************************
*
*   NUCHAR  DETERMINS AND RETURNS THE REAL VALUES OF ALL NUMBERS
*           FOUND IN 'LINE'. ALL CONNECTED SUBSTRINGS ARE ASSUMED
*           TO CONTAIN NUMBERS
*   ON ENTRY LINE    = CHARACTER STRING
*   ON EXIT  VALUE   = ARRAY OF NVALUE REAL VALUES
*
************************************************************************
      DIMENSION VALUE(40),ISTART(40)
      CHARACTER*80 LINE
      CHARACTER*1 TAB,COMMA,SPACE
      LOGICAL LEADSP
         SAVE                                                           GL0892
      DATA COMMA,SPACE/',',' '/
      TAB=CHAR(9)
*
* CLEAN OUT TABS AND COMMAS
*
      DO 10 I=1,80
   10 IF(LINE(I:I).EQ.TAB.OR.LINE(I:I).EQ.COMMA)LINE(I:I)=SPACE
*
* FIND INITIAL DIGIT OF ALL NUMBERS, CHECK FOR LEADING SPACES FOLLOWED
*     BY A CHARACTER
*
      LEADSP=.TRUE.
      NVALUE=0
      DO 20 I=1,80
         IF (LEADSP.AND.LINE(I:I).NE.SPACE) THEN
            NVALUE=NVALUE+1
            ISTART(NVALUE)=I
         ENDIF
         LEADSP=(LINE(I:I).EQ.SPACE)
   20 CONTINUE
*
* FILL NUMBER ARRAY
*
      DO 30 I=1,NVALUE
         VALUE(I)=READA(LINE,ISTART(I))
   30 CONTINUE
      RETURN
      END
      SUBROUTINE OSINV (A,N,D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION A(*)
************************************************************************
*
*    OSINV INVERTS A GENERAL SQUARE MATRIX OF ORDER UP TO MAXORB. SEE
*          DIMENSION STATEMENTS BELOW.
*
*   ON INPUT       A = GENERAL SQUARE MATRIX STORED LINEARLY.
*                  N = DIMENSION OF MATRIX A.
*                  D = VARIABLE, NOT DEFINED ON INPUT.
*
*   ON OUTPUT      A = INVERSE OF ORIGINAL A.
*                  D = DETERMINANT OF ORIGINAL A, UNLESS A WAS SINGULAR,
*                      IN WHICH CASE D = 0.0
*
************************************************************************
      DIMENSION L(MAXORB), M(MAXORB)
         SAVE                                                           GL0892
************************************************************************
*
*    IF THE VALUE OF TOL GIVEN HERE IS UNSUITABLE, IT CAN BE CHANGED.
      TOL=1.D-8
*
*
************************************************************************
      D=1.0
      NK=-N
      DO 180 K=1,N
         NK=NK+N
         L(K)=K
         M(K)=K
         KK=NK+K
         BIGA=A(KK)
         DO 20 J=K,N
            IZ=N*(J-1)
            DO 20 I=K,N
               IJ=IZ+I
C
C     10 FOLLOWS
C
               IF (ABS(BIGA)-ABS(A(IJ))) 10,20,20
   10          BIGA=A(IJ)
               L(K)=I
               M(K)=J
   20    CONTINUE
         J=L(K)
         IF (J-K) 50,50,30
   30    KI=K-N
         DO 40 I=1,N
            KI=KI+N
            HOLO=-A(KI)
            JI=KI-K+J
            A(KI)=A(JI)
   40    A(JI)=HOLO
   50    I=M(K)
         IF (I-K) 80,80,60
   60    JP=N*(I-1)
         DO 70 J=1,N
            JK=NK+J
            JI=JP+J
            HOLO=-A(JK)
            A(JK)=A(JI)
   70    A(JI)=HOLO
   80    IF (ABS(BIGA)-TOL) 90,100,100
   90    D=0.0
         RETURN
  100    DO 120 I=1,N
            IF (I-K) 110,120,110
  110       IK=NK+I
            A(IK)=A(IK)/(-BIGA)
  120    CONTINUE
         DO 150 I=1,N
            IK=NK+I
            IJ=I-N
            DO 150 J=1,N
               IJ=IJ+N
               IF (I-K) 130,150,130
  130          IF (J-K) 140,150,140
  140          KJ=IJ-I+K
               A(IJ)=A(IK)*A(KJ)+A(IJ)
  150    CONTINUE
         KJ=K-N
         DO 170 J=1,N
            KJ=KJ+N
            IF (J-K) 160,170,160
  160       A(KJ)=A(KJ)/BIGA
  170    CONTINUE
         D=D*BIGA
         A(KK)=1.0/BIGA
  180 CONTINUE
      K=N
  190 K=K-1
      IF (K) 260,260,200
  200 I=L(K)
      IF (I-K) 230,230,210
  210 JQ=N*(K-1)
      JR=N*(I-1)
      DO 220 J=1,N
         JK=JQ+J
         HOLO=A(JK)
         JI=JR+J
         A(JK)=-A(JI)
  220 A(JI)=HOLO
  230 J=M(K)
      IF (J-K) 190,190,240
  240 KI=K-N
      DO 250 I=1,N
         KI=KI+N
         HOLO=A(KI)
         JI=KI+J-K
         A(KI)=-A(JI)
  250 A(JI)=HOLO
      GO TO 190
  260 RETURN
C
      END
      SUBROUTINE PARSAV(MODE,N,M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
**********************************************************************
*
*   PARSAV SAVES AND RESTORES DATA USED IN NLLSQ GRADIENT MINIMIZATION.
*
*    IF MODE IS 0 DATA ARE RESTORED, IF 1 THEN SAVED.
*
**********************************************************************
      COMMON /DENSTY/ P(MPACK), PA(MPACK), PB(MPACK)
      COMMON /ALPARM/ ALPARM(3,MAXPAR),X0, X1, X2, ILOOP
      COMMON /ELEMTS/ ELEMNT(107)
      COMMON /KEYWRD/ KEYWRD
      COMMON /GEOSYM/ NDEP,LOCPAR(MAXPAR),IDEPFN(MAXPAR),
     1                     LOCDEP(MAXPAR)
      COMMON /TITLES/ KOMENT,TITLE
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1                NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /GEOM  / GEO(3,NUMATM)
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /NLLCOM/ DDDUM(6),EFSLST(MAXPAR),Q(MAXPAR,MAXPAR),
     1                R(MAXPAR,MAXPAR),XLAST(MAXPAR),
     2                IIIUM(7),IDUMY(2*MAXPAR*MAXPAR-19-MAXPAR*4)
      COMMON /GEOVAR/ XDUMMY(MAXPAR),NVAR,LOC(2,MAXPAR)                 IR0394
      COMMON /LOCVAR/ LOCVAR(2,MAXPAR)
      COMMON /VALVAR/ VALVAR(MAXPAR),NUMVAR
      DIMENSION IEL1(3),QQ(3), COORD(3,NUMATM)
      CHARACTER ELEMNT*2, KEYWRD*80,KOMENT*80, TITLE*80
         SAVE                                                           GL0892
      LATOM=0
      LPARAM=0
      OPEN(UNIT=9,FILE='FOR009',STATUS='UNKNOWN',FORM='UNFORMATTED')
      REWIND 9
      OPEN(UNIT=10,FILE='FOR010',STATUS='UNKNOWN',FORM='UNFORMATTED')
      REWIND 10
      IF(MODE.NE.0) GOTO 10
*
*  MODE=0: RETRIEVE DATA FROM DISK.
*
      READ(9,END=70,ERR=70)IIIUM,DDDUM,EFSLST,N,(XLAST(I),I=1,N),M
      READ(9)((Q(J,I),J=1,M),I=1,M)
      READ(9)((R(J,I),J=1,N),I=1,N)
      READ(9)(VALVAR(I),I=1,N)
      RETURN
   10 CONTINUE
      IF(MODE.EQ.1)THEN
         WRITE(6,'(//10X,'' **** TIME UP ****'')')
         WRITE(6,'(//10X,'' CURRENT VALUES OF GEOMETRIC VARIABLES'',//)'
     1)
         IF(NA(1) .EQ. 99) THEN
C
C  CONVERT FROM CARTESIAN COORDINATES TO INTERNAL
C
            DO 20 I=1,NATOMS
               DO 20 J=1,3
   20       COORD(J,I)=GEO(J,I)
            CALL XYZINT(COORD,NUMAT,NA,NB,NC,1.D0,GEO)
         ENDIF
         DEGREE=57.29577951D0
         GEO(2,1)=0.D0
         GEO(3,1)=0.D0
         GEO(1,1)=0.D0
         GEO(2,2)=0.D0
         GEO(3,2)=0.D0
         GEO(3,3)=0.D0
         IVAR=1
         NA(1)=0
         WRITE(6,'(A)')KEYWRD,KOMENT,TITLE
         DO 50 I=1,NATOMS
            DO 30 J=1,3
   30       IEL1(J)=0
   40       CONTINUE
            IF(LOC(1,IVAR).EQ.I) THEN
               IEL1(LOC(2,IVAR))=1
               IVAR=IVAR+1
               GOTO 40
            ENDIF
            IF(I.LT.4) THEN
               IEL1(3)=0
               IF(I.LT.3) THEN
                  IEL1(2)=0
                  IF(I.LT.2) THEN
                     IEL1(1)=0
                  ENDIF
               ENDIF
            ENDIF
            QQ(1)=GEO(1,I)
            QQ(2)=GEO(2,I)*DEGREE
            QQ(3)=GEO(3,I)*DEGREE
   50    WRITE(6,'(2X,A2,3(F12.6,I3),I4,2I3)')
     1    ELEMNT(LABELS(I)),(QQ(K),IEL1(K),K=1,3),NA(I),NB(I),NC(I)
         I=0
         X=0.D0
         WRITE(6,'(I4,3(F12.6,I3),I4,2I3)')
     1    I,X,I,X,I,X,I,I,I,I
         IF(NDEP.NE.0)THEN
            DO 60 I=1,NDEP
   60       WRITE(6,'(3(I4,'',''))')LOCPAR(I),IDEPFN(I),LOCDEP(I)
            WRITE(6,*)
         ENDIF
         WRITE(6,'(//10X,
     1''TO RESTART CALCULATION USE THE KEYWORD "RESTART".'')')
      ENDIF
      WRITE(9)IIIUM,DDDUM,EFSLST,N,(XLAST(I),I=1,N),M
      WRITE(9)((Q(J,I),J=1,M),I=1,M)
      WRITE(9)((R(J,I),J=1,N),I=1,N)
      WRITE(9)(VALVAR(I),I=1,N)
C*****
C     The density matrix is required by ITER upon restart .
C
      LINEAR=(NORBS*(NORBS+1))/2
      WRITE(10)(PA(I),I=1,LINEAR)
      IF(NALPHA.NE.0)WRITE(10)(PB(I),I=1,LINEAR)
C*****
      CLOSE(9)
      CLOSE(10)
      RETURN
   70 WRITE(6,'(//10X,''NO RESTART FILE EXISTS!'')')
      STOP
      END
      SUBROUTINE PATHS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
C    changed common path for portability  (IR)
      COMMON /PATHI / LATOM,LPARAM
      COMMON /PATHR / REACT(200)
c
      COMMON /GEOVAR/ XPARAM(MAXPAR), NVAR, LOC(2,MAXPAR)               IR0394
      COMMON /KEYWRD/ KEYWRD
      COMMON /TIMER / TIME0
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /ALPARM/ ALPARM(3,MAXPAR),X0, X1, X2, ILOOP
************************************************************************
*
*   PATH FOLLOWS A REACTION COORDINATE.   THE REACTION COORDINATE IS ON
*        ATOM LATOM, AND IS A DISTANCE IF LPARAM=1,
*                           AN ANGLE   IF LPARAM=2,
*                           AN DIHEDRALIF LPARAM=3.
*
************************************************************************
      DIMENSION GD(MAXPAR),XLAST(MAXPAR),MDFP(20),XDFP(20)
      CHARACTER*80 KEYWRD
      CHARACTER*10 TYPE(3)
         SAVE                                                           GL0892
      DATA TYPE / 'ANGSTROMS ','DEGREES   ','DEGREES   '/
      ILOOP=1
      IF(INDEX(KEYWRD,'RESTAR') .NE. 0) THEN
         MDFP(9)=0
         CALL DFPSAV(TOTIME,XPARAM,GD,XLAST,FUNCT1,MDFP,XDFP)
         WRITE(6,'(//10X,'' RESTARTING AT POINT'',I3)')ILOOP
      ENDIF
      IF(ILOOP.GT.1) GOTO 10
      WRITE(6,'(''  ABOUT TO ENTER FLEPO FROM PATH'')')
cmgc      TIME0=SECOND()
      CALL FLEPO(XPARAM,NVAR,FUNCT)
      WRITE(6,'(''  OPTIMIZED VALUES OF PARAMETERS, INITIAL POINT'')')
      CALL WRITMO(TIME0,FUNCT)
cmgc      TIME0=SECOND()
   10 CONTINUE
      IF(ILOOP.GT.2) GOTO 40
      GEO(LPARAM,LATOM)=REACT(2)
      IF(ILOOP.EQ.1) THEN
         X0=REACT(1)
         X1=X0
         X2=REACT(2)
         IF(X2.LT. -100.D0) STOP
         DO 20 I=1,NVAR
            ALPARM(2,I)=XPARAM(I)
   20    ALPARM(1,I)=XPARAM(I)
         ILOOP=2
      ENDIF
      CALL FLEPO(XPARAM,NVAR,FUNCT)
      RNORD=REACT(2)
      IF(LPARAM.GT.1) RNORD=RNORD*57.29577951D0
      WRITE(6,'(1X,16(''*****'')//17X,''REACTION COORDINATE = ''
     1,F12.4,2X,A10,19X//1X,16(''*****''))')RNORD,TYPE(LPARAM)
      CALL WRITMO(TIME0,FUNCT)
cmgc      TIME0=SECOND()
      DO 30 I=1,NVAR
   30 ALPARM(3,I)=XPARAM(I)
C
C   NOW FOR THE MAIN INTERPOLATION ROUTE
C
      IF(ILOOP.EQ.2)ILOOP=3
   40 CONTINUE
      LPR=ILOOP
      DO 110 ILOOP = LPR,100
C
         IF(REACT(ILOOP).LT. -100.D0) STOP
C
         RNORD=REACT(ILOOP)
         IF(LPARAM.GT.1) RNORD=RNORD*57.29577951D0
         WRITE(6,'(1X,16(''*****'')//19X,''REACTION COORDINATE = ''
     1,F12.4,2X,A10,19X//1X,16(''*****''))')RNORD,TYPE(LPARAM)
C
         X3=REACT(ILOOP)
         C3=(X0**2-X1**2)*(X1-X2)-(X1**2-X2**2)*(X0-X1)
C      WRITE(6,'(''   C3:'',F13.7)')C3
         IF (ABS(C3) .LT. 1.D-8) THEN
C
C    WE USE A LINEAR INTERPOLATION
C
            CC1=0.D0
            CC2=0.D0
         ELSE
C    WE DO A QUADRATIC INTERPOLATION
C
            CC1=(X1-X2)/C3
            CC2=(X0-X1)/C3
         ENDIF
         CB1=1.D0/(X1-X2)
         CB2=(X1**2-X2**2)*CB1
C
C    NOW TO CALCULATE THE INTERPOLATED COORDINATES
C
         DO 50 I=1,NVAR
            DELF0=ALPARM(1,I)-ALPARM(2,I)
            DELF1=ALPARM(2,I)-ALPARM(3,I)
            ACONST = CC1*DELF0-CC2*DELF1
            BCONST = CB1*DELF1-ACONST*CB2
            CCONST = ALPARM(3,I) - BCONST*X2 - ACONST*X2**2
            XPARAM(I)=CCONST+BCONST*X3+ACONST*X3**2
            ALPARM(1,I)=ALPARM(2,I)
   50    ALPARM(2,I)=ALPARM(3,I)
C
C   NOW TO CHECK THAT THE GUESSED GEOMETRY IS NOT TOO ABSURD
C
         DO 60 I=1,NVAR
   60    IF(ABS(XPARAM(I)-ALPARM(3,I)) .GT. 0.2) GOTO 70
         GOTO 90
   70    WRITE(6,'('' GEOMETRY TOO UNSTABLE FOR EXTRAPOLATION TO BE USED
     1''/ ,'' - THE LAST GEOMETRY IS BEING USED TO START THE NEXT''
     2,'' CALCULATION'')')
         DO 80 I=1,NVAR
   80    XPARAM(I)=ALPARM(3,I)
   90    CONTINUE
         X0=X1
         X1=X2
         X2=X3
         GEO(LPARAM,LATOM)=REACT(ILOOP)
         CALL FLEPO(XPARAM,NVAR,FUNCT)
         CALL WRITMO(TIME0,FUNCT)
cmgc         TIME0=SECOND()
         DO 100 I=1,NVAR
  100    ALPARM(3,I)=XPARAM(I)
  110 CONTINUE
      END
      SUBROUTINE PERM(IPERM,NELS,NMOS,MAXMOS,NPERMS)
      DIMENSION IPERM(MAXMOS,60), IADD(20), NEL(20)
         SAVE                                                           GL0892
************************************************************************
*
*  PERM PERMUTES NELS ENTITIES AMONG NMOS LOCATIONS. THE ENTITIES AND
*       LOCATIONS ARE EACH INDISTINGUISHABLE. THE PAULI EXCLUSION
*       PRINCIPLE IS FOLLOWED. THE NUMBER OF STATES PRODUCED IS GIVEN
*       BY NMOS!/(NELS!*(NMOS-NELS)!).
* ON INPUT: NELS  = NUMBER OF INDISTINGUISHABLE ENTITIES
*           NMOS  = NUMBER OF INDISTINGUISHABLE LOCATIONS
*
* ON OUTPUT IPERM = ARRAY OF PERMUTATIONS, A 0 INDICATES NO ENTITY,
*                   A 1 INDICATES AN ENTITY.
*           NPERM = NUMBER OF PERMUTATIONS.
*
************************************************************************
      IF(NELS.GT.NMOS)THEN
         WRITE(6,'('' NUMBER OF PARTICLES,'',I3,'' GREATER THAN NO. '',
     1''OF STATES,'',I3)')NELS,NMOS
         NPERMS=0
         RETURN
      ENDIF
      NPERMS=1
      DO 10 I=1,20
   10 NEL(I)=1000
      DO 20 I=1,NELS
   20 NEL(I)=1
      DO 50 I12=1-12+NELS,NMOS,NEL(12)
         IADD(12)=I12
         DO 50 I11=I12+1,NMOS,NEL(11)
            IADD(11)=I11
            DO 50 I10=I11+1,NMOS,NEL(10)
               IADD(10)=I10
               DO 50 I9=I10+1,NMOS,NEL(9)
                  IADD(9)=I9
                  DO 50 I8=I9+1,NMOS,NEL(8)
                     IADD(8)=I8
                     DO 50 I7=I8+1,NMOS,NEL(7)
                        IADD(7)=I7
                        DO 50 I6=I7+1,NMOS,NEL(6)
                           IADD(6)=I6
                           DO 50 I5=I6+1,NMOS,NEL(5)
                              IADD(5)=I5
                              DO 50 I4=I5+1,NMOS,NEL(4)
                                 IADD(4)=I4
                                 DO 50 I3=I4+1,NMOS,NEL(3)
                                    IADD(3)=I3
                                    DO 50 I2=I3+1,NMOS,NEL(2)
                                       IADD(2)=I2
                                       DO 50 I1=I2+1,NMOS,NEL(1)
                                          IADD(1)=I1
                                          DO 30 J=1,NMOS
   30                                     IPERM(J,NPERMS)=0
                                          DO 40 J=1,NELS
   40                                     IPERM(IADD(J),NPERMS)=1
                                          NPERMS=NPERMS+1
                                          IF(NPERMS.GT.61)THEN
                                             WRITE(6,'('' NUMBER OF PERM
     1UTATIONS TOO GREAT, LIMIT 60'')')
                                             GOTO 60
                                          ENDIF
   50 CONTINUE
   60 NPERMS=NPERMS-1
      RETURN
      END
      SUBROUTINE POLAR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
C**********************************************************************
C
C   POLAR SETS UP THE CALCULATION OF THE MOLECULAR ELECTRIC RESPONSE
C   PROPERTIES BY FFHPOL.
C
C**********************************************************************
      CHARACTER*2 ELEMNT
      COMMON /TITLES/ KOMENT,TITLE
      COMMON /POLVOL/ POLVOL(107)
      COMMON /KEYWRD/ KEYWRD
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1                NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /GEOVAR/ XPARAM(MAXPAR), NVAR, LOC(2,MAXPAR)               IR0394
      COMMON /TIMER / TIME0
      COMMON /ELEMTS/ ELEMNT(107)
      COMMON /CORE  / CORE(107)
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /GEOSYM/ NDEP, LOCPAR(MAXPAR), IDEPFN(MAXPAR),
     1                    LOCDEP(MAXPAR)
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /LAST  / LAST
      COMMON /COORD / COORD(3,NUMATM)
      COMMON /EULER / TVEC(3,3),IDTVEC
      COMMON /SCRACH/ RXYZ(MPACK),XDUMY(MAXPAR**2-MPACK)
      DIMENSION GRAD(MAXPAR),ROTVEC(3,3),DIPVEC(3),
     1          TEMPV(3,3)
      CHARACTER  KEYWRD*80, TYPE*7, KOMENT*80, TITLE*80
      LOGICAL LET
         SAVE                                                           GL0892
      TYPE=' MNDO  '
      LET=(INDEX(KEYWRD,'LET').NE.0)
      IF(INDEX(KEYWRD,'MINDO') .NE. 0) TYPE='MINDO/3'
      IF(INDEX(KEYWRD,'AM1') .NE. 0)    TYPE='  AM1  '
      WRITE (6,10)
   10 FORMAT ('1',20('*'),' FINITE-FIELD POLARIZABILITIES ',
     1        20('*'),//)
      CALL GMETRY(GEO,COORD)
C
C  Orient the molecule with the moments of inertia.
C  This is done to ensure a unique, reproduceable set of directions.
C  If LET is specified, the input orientation will be used.
C
      IF (.NOT.LET) THEN
         MASS = 1.0D00
         CALL AXIS (COORD,NUMAT,A,B,C,SUMW,MASS,ROTVEC)
         WRITE (6,20)
   20    FORMAT (/' ROTATION MATRIX FOR ORIENTATION OF MOLECULE:'/)
         DO 40 I = 1,3
            WRITE (6,30) (ROTVEC(I,J),J=1,3)
   30       FORMAT (5X,3F12.6)
   40    CONTINUE
C
C  ROTATE ATOMS
C
         DO 70 I = 1,NUMAT
            DO 60 J = 1,3
               SUM = 0.0D00
               DO 50 K = 1,3
                  SUM = SUM + COORD(K,I)*ROTVEC(K,J)
   50          CONTINUE
               GEO(J,I) = SUM
   60       CONTINUE
   70    CONTINUE
         DO 90 I = 1,NUMAT
            DO 80 J = 1,3
               COORD(J,I) = GEO(J,I)
   80       CONTINUE
   90    CONTINUE
         WRITE(6,'(//10X,''CARTESIAN COORDINATES '',/)')
         WRITE(6,'(4X,''NO.'',7X,''ATOM'',9X,''X'',
     1  9X,''Y'',9X,''Z'',/)')
         L=0
         DO 100 I=1,NUMAT
            IF(NAT(I).EQ.99.OR.NAT(I).EQ.107) GOTO 100
            L=L+1
            WRITE(6,'(I6,8X,A2,4X,3F10.4)')
     1           L,ELEMNT(NAT(I)),(COORD(J,L),J=1,3)
  100    CONTINUE
C
C  IF POLYMER, ROTATE TVEC
C
         IF (IDTVEC.GT.0) THEN
            DO 130 I = 1,IDTVEC
               DO 120 J = 1,3
                  SUM = 0.0D00
                  DO 110 K = 1,3
                     SUM = SUM + TVEC(K,I)*ROTVEC(K,J)
  110             CONTINUE
                  TEMPV(J,I) = SUM
  120          CONTINUE
  130       CONTINUE
            DO 150 I = 1,3
               DO 140 J = 1,IDTVEC
                  TVEC(I,J) = TEMPV(I,J)
  140          CONTINUE
  150       CONTINUE
            WRITE (6,160) ((TVEC(J,I),J=1,3),I=1,IDTVEC)
  160       FORMAT (/' NEW TRANSLATION VECTOR:'/,
     1           ' ',3(3F15.5))
         ENDIF
      ENDIF
C
      LAST=1
      NA(1)=99
C
C  SET UP THE VARIABLES IN XPARAM AND LOC, THESE ARE IN CARTESIAN
C  COORDINATES.
C
      NDEP=0
      NUMAT=0
      SUMX=0.D0
      SUMY=0.D0
      SUMZ=0.D0
      DO 180 I=1,NATOMS
         IF((LABELS(I).NE.99).AND.(LABELS(I).NE.107)) THEN
            NUMAT=NUMAT+1
            LABELS(NUMAT)=LABELS(I)
            SUMX=SUMX+COORD(1,NUMAT)
            SUMY=SUMY+COORD(2,NUMAT)
            SUMZ=SUMZ+COORD(3,NUMAT)
            DO 170 J=1,3
  170       GEO(J,NUMAT)=COORD(J,NUMAT)
         ENDIF
  180 CONTINUE
      SUMX=SUMX/NUMAT
      SUMY=SUMY/NUMAT
      SUMZ=SUMZ/NUMAT
      SUMMAX=0.D0
      ATPOL=0.D0
      DO 190 I=1,NUMAT
         IF (LABELS(I).NE.107) THEN
            ATPOL=ATPOL+POLVOL(LABELS(I))
         ENDIF
         GEO(1,I)=GEO(1,I)-SUMX
         IF(SUMMAX.LT.ABS(GEO(1,I))) SUMMAX=ABS(GEO(1,I))
         GEO(2,I)=GEO(2,I)-SUMY
         IF(SUMMAX.LT.ABS(GEO(2,I))) SUMMAX=ABS(GEO(2,I))
         GEO(3,I)=GEO(3,I)-SUMZ
         IF(SUMMAX.LT.ABS(GEO(3,I))) SUMMAX=ABS(GEO(3,I))
  190 CONTINUE
C
      NVAR=0
      NATOMS = NUMAT
      CALL COMPFG(GEO, .TRUE., HEAT0, .TRUE., GRAD, .FALSE.)
      WRITE (6,200) HEAT0
  200 FORMAT (//' ENERGY OF "REORIENTED" SYSTEM WITHOUT FIELD:',
     1        F20.10)
C
      CALL FFHPOL (HEAT0,ATPOL,DIPVEC)
C
      RETURN
      END
      SUBROUTINE QUADR(F0,F1,F2, X1,X2, A,B,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         SAVE                                                           GL0892
****************************************************
*                                                  *
*    QUADR CALCULATES THE A, B AND C IN THE EQUNS. *
*                                                  *
*     A                   =   F0                   *
*     A + B.X0 + C.X0**2  =   F1                   *
*     A + B.X2 + C.X2**2  =   F2                   *
*                                                  *
****************************************************
      C=(X2*(F1-F0)-X1*(F2-F0))/(X2*X1**2-X1*X2**2)
      B=(F1-F0-C*X1**2)/X1
      A=F0
      RETURN
      END
      DOUBLE PRECISION FUNCTION READA(STRING,ISTART)
C     FORTRAN FUNCTION TO EXTRACT NUMBER FROM STRING
C
      CHARACTER STRING*(*)
      DOUBLE PRECISION DIGIT
      LOGICAL DEFALT,EXPNNT
         SAVE                                                           GL0892
C
C     DEFINE ASCII VALUES OF NUMERIC FIELD CHARACTERS
      I0=ICHAR('0')
      I9=ICHAR('9')
      IDOT=ICHAR('.')
      INEG=ICHAR('-')
      IPOS=ICHAR('+')
      ICAPD=ICHAR('D')
      ICAPE=ICHAR('E')
      ISMLD=ICHAR('d')
      ISMLE=ICHAR('e')
C
      L=LEN(STRING)
C
C     FIND THE START OF THE NUMERIC FIELD
      DO 10 I=ISTART,L
         IADD=0
         N=ICHAR(STRING(I:I))
C
C       SIGNAL START OF NUMERIC FIELD IF DIGIT FOUND
         IF(N.GE.I0.AND.N.LE.I9)GOTO 20
C
C       ACCOUNT FOR CONSECUTIVE SIGNS [- AND(OR) +]
         IF(N.EQ.INEG.OR.N.EQ.IPOS)THEN
            IADD=IADD+1
            IF(I+IADD.GT.L)GOTO 50
            N=ICHAR(STRING(I+IADD:I+IADD))
            IF(N.GE.I0.AND.N.LE.I9)GOTO 20
         ENDIF
C
C       ACCOUNT FOR CONSECUTIVE DECIMAL POINTS (.)
         IF(N.EQ.IDOT)THEN
            IADD=IADD+1
            IF(I+IADD.GT.L)GOTO 50
            N=ICHAR(STRING(I+IADD:I+IADD))
            IF(N.GE.I0.AND.N.LE.I9)GOTO 20
         ENDIF
   10 CONTINUE
      GOTO 50
C
C     FIND THE END OF THE NUMERIC FIELD
   20 EXPNNT=.FALSE.
      DO 30 J=I+1,L
         IADD=0
         N=ICHAR(STRING(J:J))
C
C       CONTINUE SEARCH FOR END IF DIGIT FOUND
         IF(N.GE.I0.AND.N.LE.I9)GOTO 30
C
C       CONTINUE SEARCH FOR END IF SIGN FOUND AND EXPNNT TRUE
         IF(N.EQ.INEG.OR.N.EQ.IPOS)THEN
            IF(.NOT.EXPNNT)GOTO 40
            IADD=IADD+1
            IF(J+IADD.GT.L)GOTO 40
            N=ICHAR(STRING(J+IADD:J+IADD))
            IF(N.GE.I0.AND.N.LE.I9)GOTO 30
         ENDIF
         IF(N.EQ.IDOT)THEN
            IADD=IADD+1
            IF(J+IADD.GT.L)GOTO 40
            N=ICHAR(STRING(J+IADD:J+IADD))
            IF(N.GE.I0.AND.N.LE.I9)GOTO 30
            IF(N.EQ.ICAPE.OR.N.EQ.ISMLE.OR.N.EQ.ICAPD.OR.N.EQ.ISMLD)
     1    GOTO 30
         ENDIF
         IF(N.EQ.ICAPE.OR.N.EQ.ISMLE.OR.N.EQ.ICAPD.OR.N.EQ.ISMLD)THEN
            IF(EXPNNT)GOTO 40
            EXPNNT=.TRUE.
            GOTO 30
         ENDIF
         GOTO 40
   30 CONTINUE
      J=L+1
   40 N=ICHAR(STRING(J-1:J-1))
      IF(N.EQ.ICAPE.OR.N.EQ.ISMLE.OR.N.EQ.ICAPD.OR.N.EQ.ISMLD)J=J-1
C
C     FOUND THE END OF THE NUMERIC FIELD (IT RUNS 'I' THRU 'J-1')
      N=0
      N=N+INDEX(STRING(I:J-1),'e')
      N=N+INDEX(STRING(I:J-1),'E')
      N=N+INDEX(STRING(I:J-1),'d')
      N=N+INDEX(STRING(I:J-1),'D')
      IF(N.EQ.0)THEN
         READA=DIGIT(STRING(I:J-1),1)
      ELSE
         READA=DIGIT(STRING(:I+N-2),I)*1.D1**DIGIT(STRING(:J-1),I+N)
      ENDIF
      DEFALT=.FALSE.
      RETURN
C
C     DEFAULT VALUE RETURNED BECAUSE NO NUMERIC FIELD FOUND
   50 READA=0.D0
      DEFALT=.TRUE.
      RETURN
      END
      SUBROUTINE REFER
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      COMMON /REFS/ ALLREF(107,4)
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /KEYWRD/ KEYWRD
      LOGICAL ALLOK, ELEMNS(107), MIXOK, MIX
      CHARACTER*80 KEYWRD, ALLREF
         SAVE                                                           GL0892
      DATA MIX/.FALSE./
      MIXOK=(INDEX(KEYWRD,'PARASOK').NE.0)
      DO 10 I=1,102
   10 ELEMNS(I)=.FALSE.
      IF(INDEX(KEYWRD,'PM3').NE.0)THEN
         MODE=4
      ELSEIF(INDEX(KEYWRD,'AM1').NE.0)THEN
         MODE=3
      ELSEIF(INDEX(KEYWRD,'MINDO').NE.0)THEN
         MODE=2
      ELSE
         MODE=1
      ENDIF
      ALLREF(99,MODE)=' DUMMY ATOMS ARE USED; THESE DO NOT AFFECT '
     1//'THE CALCULATION'
      ALLREF(100,MODE)=' '
      DO 20 I=1,NUMAT
         J=NAT(I)
   20 ELEMNS(J)=.TRUE.
      ALLOK=.TRUE.
      DO 30 I=1,102
         IF(ELEMNS(I))THEN
            IF(I.LT.99.AND..NOT.MIX.AND.MODE.EQ.3)
     1MIX=(INDEX(ALLREF(I,3),'MNDO').NE.0)
            IF(ALLREF(I,MODE)(1:1).NE.' ')THEN
               WRITE(6,'(A,I3)')' DATA ARE NOT AVAILABLE FOR ELEMENT NO.
     1',I
               ALLOK=.FALSE.
            ELSE
               WRITE(6,'(A)')ALLREF(I,MODE)
            ENDIF
         ENDIF
   30 CONTINUE
      IF(MIX.AND..NOT.MIXOK)THEN
         WRITE(6,40)
     1  'SOME ELEMENTS HAVE BEEN SPECIFIED FOR WHICH ONLY MNDO',
     2  'PARAMETERS ARE AVAILABLE.  SUCH MIXTURES OF METHODS ARE',
     3  'VERY RISKY AND HAVE NOT BEEN FULLY TESTED.  IF YOU FEEL',
     4  'THE RISK IS WORTH WHILE - CHECK THE MANUAL FIRST - THEN',
     5  'SPECIFY "PARASOK" IN THE KEYWORDS'
         STOP
      ENDIF
      IF(ALLOK)RETURN
      WRITE(6,40)
     1 'SOME ELEMENTS HAVE BEEN SPECIFIED FOR WHICH',
     2 'NO PARAMETERS ARE AVAILABLE.  CALCULATION STOPPED.'
   40 FORMAT(/////10X,A,4(/10X,A))
      STOP
      END
      SUBROUTINE REPP(NI,NJ,RIJ,RI,CORE)
C***********************************************************************
C
C..VECTOR VERSION WRITTEN BY ERNEST R. DAVIDSON, INDIANA UNIVERSITY
C
C
C  REPP CALCULATES THE TWO-ELECTRON REPULSION INTEGRALS AND THE
C       NUCLEAR ATTRACTION INTEGRALS.
C
C     ON INPUT RIJ     = INTERATOMIC DISTANCE
C              NI      = ATOM NUMBER OF FIRST ATOM
C              NJ      = ATOM NUMBER OF SECOND ATOM
C    (REF)     ADD     = ARRAY OF GAMMA, OR TWO-ELECTRON ONE-CENTER,
C                        INTEGRALS.
C    (REF)     TORE    = ARRAY OF NUCLEAR CHARGES OF THE ELEMENTS
C    (REF)     DD      = ARRAY OF DIPOLE CHARGE SEPARATIONS
C    (REF)     QQ      = ARRAY OF QUADRUPOLE CHARGE SEPARATIONS
C
C     THE COMMON BLOCKS ARE INITIALIZED IN BLOCK-DATA, AND NEVER CHANGED
C
C    ON OUTPUT RI      = ARRAY OF TWO-ELECTRON REPULSION INTEGRALS
C              CORE    = 4 X 2 ARRAY OF ELECTRON-CORE ATTRACTION
C                        INTEGRALS
C
C
C *** THIS ROUTINE COMPUTES THE TWO-CENTRE REPULSION INTEGRALS AND THE
C *** NUCLEAR ATTRACTION INTEGRALS.
C *** THE TWO-CENTRE REPULSION INTEGRALS (OVER LOCAL COORDINATES) ARE
C *** STORED AS FOLLOWS (WHERE P-SIGMA = O,  AND P-PI = P AND P* )
C     (SS/SS)=1,   (SO/SS)=2,   (OO/SS)=3,   (PP/SS)=4,   (SS/OS)=5,
C     (SO/SO)=6,   (SP/SP)=7,   (OO/SO)=8,   (PP/SO)=9,   (PO/SP)=10,
C     (SS/OO)=11,  (SS/PP)=12,  (SO/OO)=13,  (SO/PP)=14,  (SP/OP)=15,
C     (OO/OO)=16,  (PP/OO)=17,  (OO/PP)=18,  (PP/PP)=19,  (PO/PO)=20,
C     (PP/P*P*)=21,   (P*P/P*P)=22.
C *** THE STORAGE OF THE NUCLEAR ATTRACTION INTEGRALS  CORE(KL/IJ) IS
C     (SS/)=1,   (SO/)=2,   (OO/)=3,   (PP/)=4
C     WHERE IJ=1 IF THE ORBITALS CENTRED ON ATOM I,  =2 IF ON ATOM J.
C *** NI AND NJ ARE THE ATOMIC NUMBERS OF THE TWO ELEMENTS.
C
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SI,SJ
      COMMON /MULTIP/ DD(107),QQ(107),ADD(107,3)
      COMMON /CORE/ TORE(107)
      COMMON /NATORB/ NATORB(107)
      DIMENSION RI(22),CORE(4,2)
      DIMENSION ARG(72),SQR(72)
         SAVE                                                           GL0892
      DATA  OD/1.D00/, TD/2.D00/, FD/4.D00/, ED/8.D0/
      DATA PP/0.5D00/, P2/0.25D00/, P3/0.125D00/, P4/0.0625D00/
      DATA A0/0.529167D0/ ,EV/27.21D0/, EV1/13.605D0/, EV2/6.8025D0/,
     1 EV3/3.40125D0/, EV4/1.700625D0/
C
C     ATOMIC UNITS ARE USED IN THE CALCULATION,
C     FINAL RESULTS ARE CONVERTED TO EV
C
      R=RIJ/A0
C
      SI = (NATORB(NI).GE.3)
      SJ = (NATORB(NJ).GE.3)
C
      IF ((.NOT.SI) .AND. (.NOT.SJ)) THEN
C
C     HYDROGEN - HYDROGEN  (SS/SS)
C
         AEE = PP/ADD(NI,1) + PP/ADD(NJ,1)
         AEE = AEE * AEE
         RI(1) = EV/SQRT(R*R+AEE)
         CORE(1,1 )= TORE(NJ)*RI(1)
         CORE(1,2) = TORE(NI)*RI(1)
C
      ELSE IF (SI .AND. (.NOT.SJ)) THEN
C
C     HEAVY ATOM - HYDROGEN
C
         AEE = PP/ADD(NI,1) + PP/ADD(NJ,1)
         AEE = AEE * AEE
         DA=DD(NI)
         QA=QQ(NI) * TD
         ADE = PP/ADD(NI,2) + PP/ADD(NJ,1)
         ADE = ADE * ADE
         AQE = PP/ADD(NI,3) + PP/ADD(NJ,1)
         AQE = AQE * AQE
         RSQ = R*R
         ARG(1) = RSQ + AEE
         XXX = R+DA
         ARG(2) = XXX*XXX + ADE
         XXX = R-DA
         ARG(3) = XXX*XXX + ADE
         XXX = R+QA
         ARG(4) = XXX*XXX + AQE
         XXX = R-QA
         ARG(5) = XXX*XXX + AQE
         ARG(6) = RSQ + AQE
         ARG(7) = ARG(6) + QA*QA
         DO 10 I = 1,7
            SQR(I) = SQRT(ARG(I))
   10    CONTINUE
         EE = EV/SQR(1)
         RI(1) = EE
         RI(2) = EV1/SQR(2) - EV1/SQR(3)
         RI(3) = EE + EV2/SQR(4) + EV2/SQR(5) - EV1/SQR(6)
         RI(4) = EE + EV1/SQR(7) - EV1/SQR(6)
         CORE(1,1) = TORE(NJ)*RI(1)
         CORE(1,2) = TORE(NI)*RI(1)
         CORE(2,1) = TORE(NJ)*RI(2)
         CORE(3,1) = TORE(NJ)*RI(3)
         CORE(4,1) = TORE(NJ)*RI(4)
C
      ELSE IF ((.NOT.SI).AND.SJ) THEN
C
C     HYDROGEN - HEAVY ATOM
C
         AEE = PP/ADD(NI,1) + PP/ADD(NJ,1)
         AEE = AEE * AEE
         DB=DD(NJ)
         QB=QQ(NJ) * TD
         AED = PP/ADD(NI,1) + PP/ADD(NJ,2)
         AED = AED * AED
         AEQ = PP/ADD(NI,1) + PP/ADD(NJ,3)
         AEQ = AEQ * AEQ
         RSQ = R*R
         ARG(1) = RSQ + AEE
         XXX = R-DB
         ARG(2) = XXX*XXX + AED
         XXX = R+DB
         ARG(3) = XXX*XXX + AED
         XXX = R-QB
         ARG(4) = XXX*XXX + AEQ
         XXX = R+QB
         ARG(5) = XXX*XXX + AEQ
         ARG(6) = RSQ + AEQ
         ARG(7) = ARG(6) + QB*QB
         DO 20 I = 1,7
            SQR(I) = SQRT(ARG(I))
   20    CONTINUE
         EE = EV/SQR(1)
         RI(1) = EE
         RI(5) = EV1/SQR(2)  - EV1/SQR(3)
         RI(11) = EE + EV2/SQR(4) + EV2/SQR(5) - EV1/SQR(6)
         RI(12) = EE + EV1/SQR(7) - EV1/SQR(6)
         CORE(1,1) = TORE(NJ)*RI(1)
         CORE(1,2) = TORE(NI)*RI(1)
         CORE(2,2) = TORE(NI)*RI(5)
         CORE(3,2) = TORE(NI)*RI(11)
         CORE(4,2) = TORE(NI)*RI(12)
C
      ELSE
C
C     HEAVY ATOM - HEAVY ATOM
C
C     DEFINE CHARGE SEPARATIONS.
         DA=DD(NI)
         DB=DD(NJ)
         QA=QQ(NI) * TD
         QB=QQ(NJ) * TD
C
         AEE = PP/ADD(NI,1) + PP/ADD(NJ,1)
         AEE = AEE * AEE
C
         ADE = PP/ADD(NI,2) + PP/ADD(NJ,1)
         ADE = ADE * ADE
         AQE = PP/ADD(NI,3) + PP/ADD(NJ,1)
         AQE = AQE * AQE
         AED = PP/ADD(NI,1) + PP/ADD(NJ,2)
         AED = AED * AED
         AEQ = PP/ADD(NI,1) + PP/ADD(NJ,3)
         AEQ = AEQ * AEQ
         AXX = PP/ADD(NI,2) + PP/ADD(NJ,2)
         AXX = AXX * AXX
         ADQ = PP/ADD(NI,2) + PP/ADD(NJ,3)
         ADQ = ADQ * ADQ
         AQD = PP/ADD(NI,3) + PP/ADD(NJ,2)
         AQD = AQD * AQD
         AQQ = PP/ADD(NI,3) + PP/ADD(NJ,3)
         AQQ = AQQ * AQQ
         RSQ = R * R
         ARG(1) = RSQ + AEE
         XXX = R + DA
         ARG(2) = XXX * XXX + ADE
         XXX = R - DA
         ARG(3) = XXX*XXX + ADE
         XXX = R - QA
         ARG(4) = XXX*XXX + AQE
         XXX = R + QA
         ARG(5) = XXX*XXX + AQE
         ARG(6) = RSQ + AQE
         ARG(7) = ARG(6) + QA*QA
         XXX = R-DB
         ARG(8) = XXX*XXX + AED
         XXX = R+DB
         ARG(9) = XXX*XXX + AED
         XXX = R - QB
         ARG(10) = XXX*XXX + AEQ
         XXX = R + QB
         ARG(11) = XXX*XXX + AEQ
         ARG(12) = RSQ + AEQ
         ARG(13) = ARG(12) + QB*QB
         XXX = DA-DB
         ARG(14) = RSQ + AXX + XXX*XXX
         XXX = DA+DB
         ARG(15) = RSQ + AXX + XXX*XXX
         XXX = R + DA - DB
         ARG(16) = XXX*XXX + AXX
         XXX = R - DA + DB
         ARG(17) = XXX*XXX + AXX
         XXX = R - DA - DB
         ARG(18) = XXX*XXX + AXX
         XXX = R + DA + DB
         ARG(19) = XXX*XXX + AXX
         XXX = R + DA
         ARG(20) = XXX*XXX + ADQ
         ARG(21) = ARG(20) + QB*QB
         XXX = R - DA
         ARG(22) = XXX*XXX + ADQ
         ARG(23) = ARG(22) + QB*QB
         XXX = R - DB
         ARG(24) = XXX*XXX + AQD
         ARG(25) = ARG(24) + QA*QA
         XXX = R + DB
         ARG(26) = XXX*XXX + AQD
         ARG(27) = ARG(26) + QA*QA
         XXX = R + DA - QB
         ARG(28) = XXX*XXX + ADQ
         XXX = R - DA - QB
         ARG(29) = XXX*XXX + ADQ
         XXX = R + DA + QB
         ARG(30) = XXX*XXX + ADQ
         XXX = R - DA + QB
         ARG(31) = XXX*XXX + ADQ
         XXX = R + QA - DB
         ARG(32) = XXX*XXX + AQD
         XXX = R + QA + DB
         ARG(33) = XXX*XXX + AQD
         XXX = R - QA - DB
         ARG(34) = XXX*XXX + AQD
         XXX = R - QA + DB
         ARG(35) = XXX*XXX + AQD
         ARG(36) = RSQ + AQQ
         XXX = QA - QB
         ARG(37) = ARG(36) + XXX*XXX
         XXX = QA + QB
         ARG(38) = ARG(36) + XXX*XXX
         ARG(39) = ARG(36) + QA*QA
         ARG(40) = ARG(36) + QB*QB
         ARG(41) = ARG(39) + QB*QB
         XXX = R - QB
         ARG(42) = XXX*XXX + AQQ
         ARG(43) = ARG(42) + QA*QA
         XXX = R + QB
         ARG(44) = XXX*XXX + AQQ
         ARG(45) = ARG(44) + QA*QA
         XXX = R + QA
         ARG(46) = XXX*XXX + AQQ
         ARG(47) = ARG(46) + QB*QB
         XXX = R - QA
         ARG(48) = XXX*XXX + AQQ
         ARG(49) = ARG(48) + QB*QB
         XXX = R + QA - QB
         ARG(50) = XXX*XXX + AQQ
         XXX = R + QA + QB
         ARG(51) = XXX*XXX + AQQ
         XXX = R - QA - QB
         ARG(52) = XXX*XXX + AQQ
         XXX = R - QA + QB
         ARG(53) = XXX*XXX + AQQ
         QA=QQ(NI)
         QB=QQ(NJ)
         XXX = DA - QB
         XXX = XXX*XXX
         YYY = R - QB
         YYY = YYY*YYY
         ZZZ = DA + QB
         ZZZ = ZZZ*ZZZ
         WWW = R + QB
         WWW = WWW*WWW
         ARG(54) = XXX + YYY + ADQ
         ARG(55) = XXX + WWW + ADQ
         ARG(56) = ZZZ + YYY + ADQ
         ARG(57) = ZZZ + WWW + ADQ
         XXX = QA - DB
         XXX = XXX*XXX
         YYY = QA + DB
         YYY = YYY*YYY
         ZZZ = R + QA
         ZZZ = ZZZ*ZZZ
         WWW = R - QA
         WWW = WWW*WWW
         ARG(58) = ZZZ + XXX + AQD
         ARG(59) = WWW + XXX + AQD
         ARG(60) = ZZZ + YYY + AQD
         ARG(61) = WWW + YYY + AQD
         XXX = QA - QB
         XXX = XXX*XXX
         ARG(62) = ARG(36) + TD*XXX
         YYY = QA + QB
         YYY = YYY*YYY
         ARG(63) = ARG(36) + TD*YYY
         ARG(64) = ARG(36) + TD*(QA*QA+QB*QB)
         ZZZ = R + QA - QB
         ZZZ = ZZZ*ZZZ
         ARG(65) = ZZZ + XXX + AQQ
         ARG(66) = ZZZ + YYY + AQQ
         ZZZ = R + QA + QB
         ZZZ = ZZZ*ZZZ
         ARG(67) = ZZZ + XXX + AQQ
         ARG(68) = ZZZ + YYY + AQQ
         ZZZ = R - QA - QB
         ZZZ = ZZZ*ZZZ
         ARG(69) = ZZZ + XXX + AQQ
         ARG(70) = ZZZ + YYY + AQQ
         ZZZ = R - QA + QB
         ZZZ = ZZZ*ZZZ
         ARG(71) = ZZZ + XXX + AQQ
         ARG(72) = ZZZ + YYY + AQQ
         DO 30 I = 1,72
            SQR(I) = SQRT(ARG(I))
   30    CONTINUE
         EE = EV/SQR(1)
         DZE = -EV1/SQR(2) + EV1/SQR(3)
         QZZE = EV2/SQR(4) + EV2/SQR(5) - EV1/SQR(6)
         QXXE = EV1/SQR(7) - EV1/SQR(6)
         EDZ = - EV1/SQR(8) + EV1/SQR(9)
         EQZZ  = EV2/SQR(10) + EV2/SQR(11) - EV1/SQR(12)
         EQXX  = EV1/SQR(13) - EV1/SQR(12)
         DXDX  = EV1/SQR(14) - EV1/SQR(15)
         DZDZ  = EV2/SQR(16) + EV2/SQR(17) - EV2/SQR(18) - EV2/SQR(19)
         DZQXX =  EV2/SQR(20) - EV2/SQR(21) - EV2/SQR(22) + EV2/SQR(23)
         QXXDZ =  EV2/SQR(24) - EV2/SQR(25) - EV2/SQR(26) + EV2/SQR(27)
         DZQZZ = -EV3/SQR(28) + EV3/SQR(29) - EV3/SQR(30) + EV3/SQR(31)
     1       - EV2/SQR(22) + EV2/SQR(20)
         QZZDZ = -EV3/SQR(32) + EV3/SQR(33) - EV3/SQR(34) + EV3/SQR(35)
     1       + EV2/SQR(24) - EV2/SQR(26)
         QXXQXX = EV3/SQR(37) + EV3/SQR(38) - EV2/SQR(39) - EV2/SQR(40)
     1       + EV2/SQR(36)
         QXXQYY = EV2/SQR(41) - EV2/SQR(39) - EV2/SQR(40) + EV2/SQR(36)
         QXXQZZ = EV3/SQR(43) + EV3/SQR(45) - EV3/SQR(42) - EV3/SQR(44)
     1       - EV2/SQR(39) + EV2/SQR(36)
         QZZQXX = EV3/SQR(47) + EV3/SQR(49) - EV3/SQR(46) - EV3/SQR(48)
     1       - EV2/SQR(40) + EV2/SQR(36)
         QZZQZZ = EV4/SQR(50) + EV4/SQR(51) + EV4/SQR(52) + EV4/SQR(53)
     1       - EV3/SQR(48) - EV3/SQR(46) - EV3/SQR(42) - EV3/SQR(44)
     2       + EV2/SQR(36)
         DXQXZ = -EV2/SQR(54) + EV2/SQR(55) + EV2/SQR(56) - EV2/SQR(57)
         QXZDX = -EV2/SQR(58) + EV2/SQR(59) + EV2/SQR(60) - EV2/SQR(61)
         QXYQXY = EV2/SQR(62) + EV2/SQR(63) - EV1/SQR(64)
         QXZQXZ = EV3/SQR(65) - EV3/SQR(67) - EV3/SQR(69) + EV3/SQR(71)
     1       - EV3/SQR(66) + EV3/SQR(68) + EV3/SQR(70) - EV3/SQR(72)
         RI(1) = EE
         RI(2) = -DZE
         RI(3) = EE + QZZE
         RI(4) = EE + QXXE
         RI(5) = -EDZ
         RI(6) = DZDZ
         RI(7) = DXDX
         RI(8) = -EDZ -QZZDZ
         RI(9) = -EDZ -QXXDZ
         RI(10) = -QXZDX
         RI(11) =  EE + EQZZ
         RI(12) =  EE + EQXX
         RI(13) = -DZE -DZQZZ
         RI(14) = -DZE -DZQXX
         RI(15) = -DXQXZ
         RI(16) = EE +EQZZ +QZZE +QZZQZZ
         RI(17) = EE +EQZZ +QXXE +QXXQZZ
         RI(18) = EE +EQXX +QZZE +QZZQXX
         RI(19) = EE +EQXX +QXXE +QXXQXX
         RI(20) = QXZQXZ
         RI(21) = EE +EQXX +QXXE +QXXQYY
         RI(22) = PP * (QXXQXX -QXXQYY)
C
C     CALCULATE CORE-ELECTRON ATTRACTIONS.
C
         CORE(1,1) = TORE(NJ)*RI(1)
         CORE(2,1) = TORE(NJ)*RI(2)
         CORE(3,1) = TORE(NJ)*RI(3)
         CORE(4,1) = TORE(NJ)*RI(4)
         CORE(1,2) = TORE(NI)*RI(1)
         CORE(2,2) = TORE(NI)*RI(5)
         CORE(3,2) = TORE(NI)*RI(11)
         CORE(4,2) = TORE(NI)*RI(12)
C
      END IF
C
      RETURN
C
      END
      SUBROUTINE ROTAT(COORD,I,J,IX,RIJ,DEL1,IDX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /DERIVS/ DS(16),DG(22),DR(100),TDX(3),TDY(3),TDZ(3)
      COMMON /EXTRA/  G(22),TX(3),TY(3),TZ(3)
      DIMENSION COORD(3,25)
         SAVE                                                           GL0892
      XD=COORD(1,I)-COORD(1,J)
      YD=COORD(2,I)-COORD(2,J)
      ZD=COORD(3,I)-COORD(3,J)
      RXY=SQRT(XD*XD+YD*YD)
      RYZ=SQRT(YD*YD+ZD*ZD)
      RZX=SQRT(ZD*ZD+XD*XD)
      DO 10 IJK=1,3
         TX(IJK)=0.0D0
         TY(IJK)=0.0D0
         TZ(IJK)=0.0D0
         TDX(IJK)=0.0D0
         TDY(IJK)=0.0D0
         TDZ(IJK)=0.0D0
   10 CONTINUE
      IF(RXY.LT.1.0D-4) THEN
C   MOLECULAR Z AXIS IS PARALLEL TO DIATOMIC Z AXIS
         TX(3)=1.0D0
         IF(ZD.LT.0.0D0) TX(3)=-1.0D0
         TY(2)=1.0D0
         TZ(1)=TX(3)
         IF(IDX.EQ.1) RETURN
         IF(IX.EQ.1) TDX(1)=1.0D0/RIJ
         IF(IX.EQ.2) TDX(2)=1.0D0/RIJ
         IF(IX.EQ.1) TDZ(3)=-1.0D0/RIJ
         IF(IX.EQ.2) TDY(3)=-TX(3)/RIJ
      ELSEIF(RYZ.LT.1.0D-4) THEN
C   MOLECULAR X AXIS IS PARALLEL TO DIATOMIC Z AXIS
         TX(1)=1.0D0
         IF(XD.LT.0.0D0) TX(1)=-1.0D0
         TY(2)=TX(1)
         TZ(3)=1.0D0
         IF(IDX.EQ.1) RETURN
         IF(IX.EQ.2) TDX(2)=1.0D0/RIJ
         IF(IX.EQ.3) TDX(3)=1.0D0/RIJ
         IF(IX.EQ.2) TDY(1)=-1.0D0/RIJ
         IF(IX.EQ.3) TDZ(1)=-TX(1)/RIJ
      ELSEIF(RZX.LT.1.0D-4) THEN
C   MOLECULAR Y AXIS IS PARALLEL TO DIATOMIC Z AXIS
         TX(2)=1.0D0
         IF(YD.LT.0.0D0) TX(2)=-1.0D0
         TY(1)=-TX(2)
         TZ(3)=1.0D0
         IF(IDX.EQ.1) RETURN
         IF(IX.EQ.1) TDX(1)=1.0D0/RIJ
         IF(IX.EQ.3) TDX(3)=1.0D0/RIJ
         IF(IX.EQ.1) TDY(2)=1.0D0/RIJ
         IF(IX.EQ.3) TDZ(2)=-TX(2)/RIJ
      ELSE
         TX(1)=XD/RIJ
         TX(2)=YD/RIJ
         TX(3)=ZD/RIJ
         TZ(3)=RXY/RIJ
         TY(1)=-TX(2)*SIGN(+1.0D0,TX(1))/TZ(3)
         TY(2)=ABS(TX(1)/TZ(3))
         TY(3)=0.0D0
         TZ(1)=-TX(1)*TX(3)/TZ(3)
         TZ(2)=-TX(2)*TX(3)/TZ(3)
         IF(IDX.EQ.1) RETURN
         TERM=DEL1/(RIJ*RIJ)
         IF(IX.EQ.1)THEN
            TDX(1)=1.0D0/RIJ-TX(1)*TERM
            TDX(2)=-TX(2)*TERM
            TDX(3)=-TX(3)*TERM
            TDZ(3)=TX(1)/RXY-TZ(3)*TERM
         ELSEIF(IX.EQ.2) THEN
            TDX(1)=-TX(1)*TERM
            TDX(2)=1.0D0/RIJ-TX(2)*TERM
            TDX(3)=-TX(3)*TERM
            TDZ(3)=TX(2)/RXY-TZ(3)*TERM
         ELSEIF(IX.EQ.3)THEN
            TDX(1)=-TX(1)*TERM
            TDX(2)=-TX(2)*TERM
            TDX(3)=1.0D0/RIJ-TX(3)*TERM
            TDZ(3)=-TZ(3)*TERM
         ENDIF
         TDY(1)=-TDX(2)/TZ(3)+TX(2)*TDZ(3)/TZ(3)**2
         IF(TX(1).LT.0.0D0) TDY(1)=-TDY(1)
         TDY(2)=TDX(1)/TZ(3)-TX(1)*TDZ(3)/TZ(3)**2
         IF(TX(1).LT.0.0D0) TDY(2)=-TDY(2)
         TDY(3)=0.0D0
         TDZ(1)=-TX(3)*TDX(1)/TZ(3)-TX(1)*TDX(3)/TZ(3)
     1 +TX(1)*TX(3)*TDZ(3)/TZ(3)**2
         TDZ(2)=-TX(3)*TDX(2)/TZ(3)-TX(2)*TDX(3)/TZ(3)
     1 +TX(2)*TX(3)*TDZ(3)/TZ(3)**2
      ENDIF
      RETURN
      END
      SUBROUTINE ROTATE (NI,NJ,XI,XJ,W,KR,E1B,E2A,ENUC,CUTOFF)
C***********************************************************************
C
C..IMPROVED SCALAR VERSION
C..WRITTEN BY ERNEST R. DAVIDSON, INDIANA UNIVERSITY.
C
C
C   ROTATE CALCULATES THE TWO-PARTICLE INTERACTIONS.
C
C   ON INPUT  NI     = ATOMIC NUMBER OF FIRST ATOM.
C             NJ     = ATOMIC NUMBER OF SECOND ATOM.
C             XI     = COORDINATE OF FIRST ATOM.
C             XJ     = COORDINATE OF SECOND ATOM.
C
C ON OUTPUT W      = ARRAY OF TWO-ELECTRON REPULSION INTEGRALS.
C           E1B,E2A= ARRAY OF ELECTRON-NUCLEAR ATTRACTION INTEGRALS,
C                    E1B = ELECTRON ON ATOM NI ATTRACTING NUCLEUS OF NJ.
C           ENUC   = NUCLEAR-NUCLEAR REPULSION TERM.
C
C
C *** THIS ROUTINE COMPUTES THE REPULSION AND NUCLEAR ATTRACTION
C     INTEGRALS OVER MOLECULAR-FRAME COORDINATES.  THE INTEGRALS OVER
C     LOCAL FRAME COORDINATES ARE EVALUATED BY SUBROUTINE REPP AND
C     STORED AS FOLLOWS (WHERE P-SIGMA = O,   AND P-PI = P AND P* )
C     IN RI
C     (SS/SS)=1,   (SO/SS)=2,   (OO/SS)=3,   (PP/SS)=4,   (SS/OS)=5,
C     (SO/SO)=6,   (SP/SP)=7,   (OO/SO)=8,   (PP/SO)=9,   (PO/SP)=10,
C     (SS/OO)=11,  (SS/PP)=12,  (SO/OO)=13,  (SO/PP)=14,  (SP/OP)=15,
C     (OO/OO)=16,  (PP/OO)=17,  (OO/PP)=18,  (PP/PP)=19,  (PO/PO)=20,
C     (PP/P*P*)=21,   (P*P/P*P)=22.
C
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*80 KEYWRD
      LOGICAL SI,SJ, AM1PM3, ANALYT, AM1
      COMMON /NATORB/ NATORB(107)
      COMMON /TWOEL3/ F03(107)
      COMMON /ALPHA3/ ALP3(153)
      COMMON /ALPHA / ALP(107)
      COMMON /CORE  / TORE(107)
      COMMON /IDEAS / FN1(107,10),FN2(107,10),FN3(107,10)
      COMMON /ALPTM / ALPTM(30), EMUDTM(30)
      COMMON /ROTDUM/ CSS1,CSP1,CPPS1,CPPP1,CSS2,CSP2,CPPS2,CPPP2
      COMMON /ROTDU2/ X(3),Y(3),Z(3)
      COMMON /KEYWRD/ KEYWRD
      DIMENSION XI(3),XJ(3),W(100),E1B(10),E2A(10)
      DIMENSION RI(22),CORE(4,2), BORON1(3,4), BORON2(3,4), BORON3(3,4)
         SAVE                                                           GL0892
      EQUIVALENCE (CORE(1,1),CSS1)

      DATA BORON1/  0.182613D0,  0.118587D0, -0.073280D0,
     1              0.412253D0, -0.149917D0,  0.000000D0,
     2              0.261751D0,  0.050275D0,  0.000000D0,
     3              0.359244D0,  0.074729D0,  0.000000D0/
      DATA BORON2/  6.D0,  6.D0,  5.D0,
     1             10.D0,  6.D0,  0.D0,
     2              8.D0,  5.D0,  0.D0,
     3              9.D0,  9.D0,  0.D0/
      DATA BORON3/  0.727592D0,  1.466639D0,  1.570975D0,
     1              0.832586D0,  1.186220D0,  0.000000D0,
     2              1.063995D0,  1.936492D0,  0.000000D0,
     3              0.819351D0,  1.574414D0,  0.000000D0/
C
C INITIALIZE ITYPE - PROGRAM ASSUMES ITYPE IS PRESERVED BETWEEN CALLS
      DATA ITYPE /0/
C
      IF (ITYPE.EQ.0) THEN
         IF(INDEX(KEYWRD,'MINDO') .NE. 0) THEN
            ITYPE=1
         ELSE
            AM1   = (INDEX(KEYWRD,'AM1') .NE. 0)
            AM1PM3= (INDEX(KEYWRD,'PM3')+INDEX(KEYWRD,'AM1') .NE. 0)
            ITYPE=2
         ENDIF
         ANALYT=(INDEX(KEYWRD,'ANALYT') .NE. 0)
         IF(ANALYT)THEN
            OPEN(UNIT=2,STATUS='SCRATCH',FORM='UNFORMATTED')
            REWIND 2
         ENDIF
      ENDIF
C
      X(1)=XI(1)-XJ(1)
      X(2)=XI(2)-XJ(2)
      X(3)=XI(3)-XJ(3)
      RIJ=X(1)*X(1)+X(2)*X(2)+X(3)*X(3)
      IF (RIJ.LT.0.002D0) THEN
C
C     SMALL RIJ CASE
C
         DO 10 I=1,10
            E1B(I)=0.D0
            E2A(I)=0.D0
   10    CONTINUE
         W(KR)=0.D0
         ENUC=0.D0
C
      ELSE IF (ITYPE.EQ.1) THEN
C
C     MINDO CASE
C
         SUM=14.399D0/SQRT(RIJ+(7.1995D0/F03(NI)+7.1995D0/F03(NJ))**2)
         W(1)=SUM
         KR=KR+1
         DO 20 L=1,10
            E1B(L)=0.D0
            E2A(L)=0.D0
   20    CONTINUE
         E1B(1) = -SUM*TORE(NJ)
         E1B(3) = E1B(1)
         E1B(6) = E1B(1)
         E1B(10)= E1B(1)
         E2A(1) = -SUM*TORE(NI)
         E2A(3) = E2A(1)
         E2A(6) = E2A(1)
         E2A(10)= E2A(1)
         II = MAX(NI,NJ)
         NBOND = (II*(II-1))/2+NI+NJ-II
         RIJ = SQRT(RIJ)
         IF(NBOND.LT.154) THEN
            IF(NBOND.EQ.22 .OR. NBOND .EQ. 29) THEN
C              NBOND = 22 IS C-H CASE
C              NBOND = 29 IS N-H CASE
               SCALE=ALP3(NBOND)*EXP(-RIJ)
            ELSE
C              NBOND < 154  IS NI < 18 AND NJ < 18 CASE
               SCALE=EXP(-ALP3(NBOND)*RIJ)
            ENDIF
         ELSE
C              NBOND > 154 INVOLVES NI OR NJ > 18
            SCALE = 0
            IF(NATORB(NI).EQ.0) SCALE=      EXP(-ALP(NI)*RIJ)
            IF(NATORB(NJ).EQ.0) SCALE=SCALE+EXP(-ALP(NI)*RIJ)
         ENDIF
         IF (ABS(TORE(NI)).GT.20.D0 .AND. ABS(TORE(NJ)).GT.20.D0) THEN
            ENUC=0.D0
         ELSE IF (RIJ.LT.1.D0 .AND. NATORB(NI)*NATORB(NJ).EQ.0) THEN
            ENUC=0.D0
         ELSE
            ENUC = TORE(NI)*TORE(NJ)*SUM
     1       + ABS(TORE(NI)*TORE(NJ)*(14.399D0/RIJ-SUM)*SCALE)
         ENDIF
C
C     MNDO AND AM1 CASES
C
C *** THE REPULSION INTEGRALS OVER MOLECULAR FRAME (W) ARE STORED IN THE
C     ORDER IN WHICH THEY WILL LATER BE USED.  IE.  (I,J/K,L) WHERE
C     J.LE.I  AND  L.LE.K     AND L VARIES MOST RAPIDLY AND I LEAST
C     RAPIDLY.  (ANTI-NORMAL COMPUTER STORAGE)
C
      ELSE
C
         RIJX = SQRT(RIJ)
         RIJ = MIN(RIJX,CUTOFF)
C
C *** COMPUTE INTEGRALS IN DIATOMIC FRAME
C
         CALL REPP(NI,NJ,RIJ,RI,CORE)
         IF(ANALYT)WRITE(2)(RI(I),I=1,22)
C
         GAM = RI(1)
         A=1.D0/RIJX
         X(1) = X(1)*A
         X(2) = X(2)*A
         X(3) = X(3)*A
         IF (ABS(X(3)).GT.0.99999999D0) THEN
            X(3) = SIGN(1.D0,X(3))
            Y(1) = 0.D0
            Y(2) = 1.D0
            Y(3) = 0.D0
            Z(1) = 1.D0
            Z(2) = 0.D0
            Z(3) = 0.D0
         ELSE
            Z(3)=SQRT(1.D0-X(3)*X(3))
            A=1.D0/Z(3)
            Y(1)=-A*X(2)*SIGN(1.D0,X(1))
            Y(2)=ABS(A*X(1))
            Y(3)=0.D0
            Z(1)=-A*X(1)*X(3)
            Z(2)=-A*X(2)*X(3)
         ENDIF
         SI = (NATORB(NI).GT.1)
         SJ = (NATORB(NJ).GT.1)
         IF ( SI .OR. SJ) THEN
            XX11 = X(1)*X(1)
            XX21 = X(2)*X(1)
            XX22 = X(2)*X(2)
            XX31 = X(3)*X(1)
            XX32 = X(3)*X(2)
            XX33 = X(3)*X(3)
            YY11 = Y(1)*Y(1)
            YY21 = Y(2)*Y(1)
            YY22 = Y(2)*Y(2)
            ZZ11 = Z(1)*Z(1)
            ZZ21 = Z(2)*Z(1)
            ZZ22 = Z(2)*Z(2)
            ZZ31 = Z(3)*Z(1)
            ZZ32 = Z(3)*Z(2)
            ZZ33 = Z(3)*Z(3)
            YYZZ11 = YY11+ZZ11
            YYZZ21 = YY21+ZZ21
            YYZZ22 = YY22+ZZ22
            XY11 = 2.D0*X(1)*Y(1)
            XY21 =      X(1)*Y(2)+X(2)*Y(1)
            XY22 = 2.D0*X(2)*Y(2)
            XY31 =      X(3)*Y(1)
            XY32 =      X(3)*Y(2)
            XZ11 = 2.D0*X(1)*Z(1)
            XZ21 =      X(1)*Z(2)+X(2)*Z(1)
            XZ22 = 2.D0*X(2)*Z(2)
            XZ31 =      X(1)*Z(3)+X(3)*Z(1)
            XZ32 =      X(2)*Z(3)+X(3)*Z(2)
            XZ33 = 2.D0*X(3)*Z(3)
            YZ11 = 2.D0*Y(1)*Z(1)
            YZ21 =      Y(1)*Z(2)+Y(2)*Z(1)
            YZ22 = 2.D0*Y(2)*Z(2)
            YZ31 =      Y(1)*Z(3)
            YZ32 =      Y(2)*Z(3)
         ENDIF
C
C     (S S/S S)
         W(1)=RI(1)
         KI = 1
         IF (SJ) THEN
C     (S S/PX S)
            W(2)=RI(5)*X(1)
C     (S S/PX PX)
            W(3)=RI(11)*XX11+RI(12)*YYZZ11
C     (S S/PY S)
            W(4)=RI(5)*X(2)
C     (S S/PY PX)
            W(5)=RI(11)*XX21+RI(12)*YYZZ21
C     (S S/PY PY)
            W(6)=RI(11)*XX22+RI(12)*YYZZ22
C     (S S/PZ S)
            W(7)=RI(5)*X(3)
C     (S S/PZ PX)
            W(8)=RI(11)*XX31+RI(12)*ZZ31
C     (S S/PZ PY)
            W(9)=RI(11)*XX32+RI(12)*ZZ32
C     (S S/PZ PZ)
            W(10)=RI(11)*XX33+RI(12)*ZZ33
            KI = 10
         ENDIF
C
         IF (SI) THEN
C     (PX S/S S)
            W(11)=RI(2)*X(1)
            IF (SJ) THEN
C     (PX S/PX S)
               W(12)=RI(6)*XX11+RI(7)*YYZZ11
C     (PX S/PX PX)
               W(13)=X(1)*(RI(13)*XX11+RI(14)*YYZZ11)
     1           +RI(15)*(Y(1)*XY11+Z(1)*XZ11)
C     (PX S/PY S)
               W(14)=RI(6)*XX21+RI(7)*YYZZ21
C     (PX S/PY PX)
               W(15)=X(1)*(RI(13)*XX21+RI(14)*YYZZ21)
     1           +RI(15)*(Y(1)*XY21+Z(1)*XZ21)
C     (PX S/PY PY)
               W(16)=X(1)*(RI(13)*XX22+RI(14)*YYZZ22)
     1           +RI(15)*(Y(1)*XY22+Z(1)*XZ22)
C     (PX S/PZ S)
               W(17)=RI(6)*XX31+RI(7)*ZZ31
C     (PX S/PZ PX)
               W(18)=X(1)*(RI(13)*XX31+RI(14)*ZZ31)
     1           +RI(15)*(Y(1)*XY31+Z(1)*XZ31)
C     (PX S/PZ PY)
               W(19)=X(1)*(RI(13)*XX32+RI(14)*ZZ32)
     1           +RI(15)*(Y(1)*XY32+Z(1)*XZ32)
C     (PX S/PZ PZ)
               W(20)=X(1)*(RI(13)*XX33+RI(14)*ZZ33)
     1           +RI(15)*(          Z(1)*XZ33)
C     (PX PX/S S)
               W(21)=RI(3)*XX11+RI(4)*YYZZ11
C     (PX PX/PX S)
               W(22)=X(1)*(RI(8)*XX11+RI(9)*YYZZ11)
     1           +RI(10)*(Y(1)*XY11+Z(1)*XZ11)
C     (PX PX/PX PX)
               W(23) =
     1     (RI(16)*XX11+RI(17)*YYZZ11)*XX11+RI(18)*XX11*YYZZ11
     2     +RI(19)*(YY11*YY11+ZZ11*ZZ11)
     3     +RI(20)*(XY11*XY11+XZ11*XZ11)
     4     +RI(21)*(YY11*ZZ11+ZZ11*YY11)
     5     +RI(22)*YZ11*YZ11
C     (PX PX/PY S)
               W(24)=X(2)*(RI(8)*XX11+RI(9)*YYZZ11)
     1           +RI(10)*(Y(2)*XY11+Z(2)*XZ11)
C     (PX PX/PY PX)
               W(25) =
     1     (RI(16)*XX11+RI(17)*YYZZ11)*XX21+RI(18)*XX11*YYZZ21
     2     +RI(19)*(YY11*YY21+ZZ11*ZZ21)
     3     +RI(20)*(XY11*XY21+XZ11*XZ21)
     4     +RI(21)*(YY11*ZZ21+ZZ11*YY21)
     5     +RI(22)*YZ11*YZ21
C     (PX PX/PY PY)
               W(26) =
     1     (RI(16)*XX11+RI(17)*YYZZ11)*XX22+RI(18)*XX11*YYZZ22
     2     +RI(19)*(YY11*YY22+ZZ11*ZZ22)
     3     +RI(20)*(XY11*XY22+XZ11*XZ22)
     4     +RI(21)*(YY11*ZZ22+ZZ11*YY22)
     5     +RI(22)*YZ11*YZ22
C     (PX PX/PZ S)
               W(27)=X(3)*(RI(8)*XX11+RI(9)*YYZZ11)
     1           +RI(10)*(         +Z(3)*XZ11)
C     (PX PX/PZ PX)
               W(28) =
     1      (RI(16)*XX11+RI(17)*YYZZ11)*XX31
     2     +(RI(18)*XX11+RI(19)*ZZ11+RI(21)*YY11)*ZZ31
     3     +RI(20)*(XY11*XY31+XZ11*XZ31)
     4     +RI(22)*YZ11*YZ31
C     (PX PX/PZ PY)
               W(29) =
     1      (RI(16)*XX11+RI(17)*YYZZ11)*XX32
     2     +(RI(18)*XX11+RI(19)*ZZ11+RI(21)*YY11)*ZZ32
     3     +RI(20)*(XY11*XY32+XZ11*XZ32)
     4     +RI(22)*YZ11*YZ32
C     (PX PX/PZ PZ)
               W(30) =
     1      (RI(16)*XX11+RI(17)*YYZZ11)*XX33
     2     +(RI(18)*XX11+RI(19)*ZZ11+RI(21)*YY11)*ZZ33
     3     +RI(20)*XZ11*XZ33
C     (PY S/S S)
               W(31)=RI(2)*X(2)
C     (PY S/PX S)
               W(32)=RI(6)*XX21+RI(7)*YYZZ21
C     (PY S/PX PX)
               W(33)=X(2)*(RI(13)*XX11+RI(14)*YYZZ11)
     1           +RI(15)*(Y(2)*XY11+Z(2)*XZ11)
C     (PY S/PY S)
               W(34)=RI(6)*XX22+RI(7)*YYZZ22
C     (PY S/PY PX)
               W(35)=X(2)*(RI(13)*XX21+RI(14)*YYZZ21)
     1           +RI(15)*(Y(2)*XY21+Z(2)*XZ21)
C     (PY S/PY PY)
               W(36)=X(2)*(RI(13)*XX22+RI(14)*YYZZ22)
     1           +RI(15)*(Y(2)*XY22+Z(2)*XZ22)
C     (PY S/PZ S)
               W(37)=RI(6)*XX32+RI(7)*ZZ32
C     (PY S/PZ PX)
               W(38)=X(2)*(RI(13)*XX31+RI(14)*ZZ31)
     1           +RI(15)*(Y(2)*XY31+Z(2)*XZ31)
C     (PY S/PZ PY)
               W(39)=X(2)*(RI(13)*XX32+RI(14)*ZZ32)
     1           +RI(15)*(Y(2)*XY32+Z(2)*XZ32)
C     (PY S/PZ PZ)
               W(40)=X(2)*(RI(13)*XX33+RI(14)*ZZ33)
     1           +RI(15)*(         +Z(2)*XZ33)
C     (PY PX/S S)
               W(41)=RI(3)*XX21+RI(4)*YYZZ21
C     (PY PX/PX S)
               W(42)=X(1)*(RI(8)*XX21+RI(9)*YYZZ21)
     1           +RI(10)*(Y(1)*XY21+Z(1)*XZ21)
C     (PY PX/PX PX)
               W(43) =
     1     (RI(16)*XX21+RI(17)*YYZZ21)*XX11+RI(18)*XX21*YYZZ11
     2     +RI(19)*(YY21*YY11+ZZ21*ZZ11)
     3     +RI(20)*(XY21*XY11+XZ21*XZ11)
     4     +RI(21)*(YY21*ZZ11+ZZ21*YY11)
     5     +RI(22)*YZ21*YZ11
C     (PY PX/PY S)
               W(44)=X(2)*(RI(8)*XX21+RI(9)*YYZZ21)
     1           +RI(10)*(Y(2)*XY21+Z(2)*XZ21)
C     (PY PX/PY PX)
               W(45) =
     1     (RI(16)*XX21+RI(17)*YYZZ21)*XX21+RI(18)*XX21*YYZZ21
     2     +RI(19)*(YY21*YY21+ZZ21*ZZ21)
     3     +RI(20)*(XY21*XY21+XZ21*XZ21)
     4     +RI(21)*(YY21*ZZ21+ZZ21*YY21)
     5     +RI(22)*YZ21*YZ21
C     (PY PX/PY PY)
               W(46) =
     1     (RI(16)*XX21+RI(17)*YYZZ21)*XX22+RI(18)*XX21*YYZZ22
     2     +RI(19)*(YY21*YY22+ZZ21*ZZ22)
     3     +RI(20)*(XY21*XY22+XZ21*XZ22)
     4     +RI(21)*(YY21*ZZ22+ZZ21*YY22)
     5     +RI(22)*YZ21*YZ22
C     (PY PX/PZ S)
               W(47)=X(3)*(RI(8)*XX21+RI(9)*YYZZ21)
     1           +RI(10)*(         +Z(3)*XZ21)
C      (PY PX/PZ PX)
               W(48) =
     1     (RI(16)*XX21+RI(17)*YYZZ21)*XX31
     2     +(RI(18)*XX21+RI(19)*ZZ21+RI(21)*YY21)*ZZ31
     3     +RI(20)*(XY21*XY31+XZ21*XZ31)
     4     +RI(22)*YZ21*YZ31
C      (PY PX/PZ PY)
               W(49) =
     1     (RI(16)*XX21+RI(17)*YYZZ21)*XX32
     2     +(RI(18)*XX21+RI(19)*ZZ21+RI(21)*YY21)*ZZ32
     3     +RI(20)*(XY21*XY32+XZ21*XZ32)
     4     +RI(22)*YZ21*YZ32
C      (PY PX/PZ PZ)
               W(50) =
     1     (RI(16)*XX21+RI(17)*YYZZ21)*XX33
     2     +(RI(18)*XX21+RI(19)*ZZ21+RI(21)*YY21)*ZZ33
     3     +RI(20)*XZ21*XZ33
C     (PY PY/S S)
               W(51)=RI(3)*XX22+RI(4)*YYZZ22
C     (PY PY/PX S)
               W(52)=X(1)*(RI(8)*XX22+RI(9)*YYZZ22)
     1           +RI(10)*(Y(1)*XY22+Z(1)*XZ22)
C      (PY PY/PX PX)
               W(53) =
     1     (RI(16)*XX22+RI(17)*YYZZ22)*XX11+RI(18)*XX22*YYZZ11
     2     +RI(19)*(YY22*YY11+ZZ22*ZZ11)
     3     +RI(20)*(XY22*XY11+XZ22*XZ11)
     4     +RI(21)*(YY22*ZZ11+ZZ22*YY11)
     5     +RI(22)*YZ22*YZ11
C     (PY PY/PY S)
               W(54)=X(2)*(RI(8)*XX22+RI(9)*YYZZ22)
     1           +RI(10)*(Y(2)*XY22+Z(2)*XZ22)
C      (PY PY/PY PX)
               W(55) =
     1     (RI(16)*XX22+RI(17)*YYZZ22)*XX21+RI(18)*XX22*YYZZ21
     2     +RI(19)*(YY22*YY21+ZZ22*ZZ21)
     3     +RI(20)*(XY22*XY21+XZ22*XZ21)
     4     +RI(21)*(YY22*ZZ21+ZZ22*YY21)
     5     +RI(22)*YZ22*YZ21
C      (PY PY/PY PY)
               W(56) =
     1     (RI(16)*XX22+RI(17)*YYZZ22)*XX22+RI(18)*XX22*YYZZ22
     2     +RI(19)*(YY22*YY22+ZZ22*ZZ22)
     3     +RI(20)*(XY22*XY22+XZ22*XZ22)
     4     +RI(21)*(YY22*ZZ22+ZZ22*YY22)
     5     +RI(22)*YZ22*YZ22
C     (PY PY/PZ S)
               W(57)=X(3)*(RI(8)*XX22+RI(9)*YYZZ22)
     1           +RI(10)*(         +Z(3)*XZ22)
C      (PY PY/PZ PX)
               W(58) =
     1     (RI(16)*XX22+RI(17)*YYZZ22)*XX31
     2     +(RI(18)*XX22+RI(19)*ZZ22+RI(21)*YY22)*ZZ31
     3     +RI(20)*(XY22*XY31+XZ22*XZ31)
     4     +RI(22)*YZ22*YZ31
C      (PY PY/PZ PY)
               W(59) =
     1     (RI(16)*XX22+RI(17)*YYZZ22)*XX32
     2     +(RI(18)*XX22+RI(19)*ZZ22+RI(21)*YY22)*ZZ32
     3     +RI(20)*(XY22*XY32+XZ22*XZ32)
     4     +RI(22)*YZ22*YZ32
C      (PY PY/PZ PZ)
               W(60) =
     1     (RI(16)*XX22+RI(17)*YYZZ22)*XX33
     2     +(RI(18)*XX22+RI(19)*ZZ22+RI(21)*YY22)*ZZ33
     3     +RI(20)*XZ22*XZ33
C     (PZ S/SS)
               W(61)=RI(2)*X(3)
C     (PZ S/PX S)
               W(62)=RI(6)*XX31+RI(7)*ZZ31
C     (PZ S/PX PX)
               W(63)=X(3)*(RI(13)*XX11+RI(14)*YYZZ11)
     1           +RI(15)*(         +Z(3)*XZ11)
C     (PZ S/PY S)
               W(64)=RI(6)*XX32+RI(7)*ZZ32
C     (PZ S/PY PX)
               W(65)=X(3)*(RI(13)*XX21+RI(14)*YYZZ21)
     1           +RI(15)*(         +Z(3)*XZ21)
C     (PZ S/PY PY)
               W(66)=X(3)*(RI(13)*XX22+RI(14)*YYZZ22)
     1           +RI(15)*(         +Z(3)*XZ22)
C     (PZ S/PZ S)
               W(67)=RI(6)*XX33+RI(7)*ZZ33
C     (PZ S/PZ PX)
               W(68)=X(3)*(RI(13)*XX31+RI(14)*ZZ31)
     1           +RI(15)*(         +Z(3)*XZ31)
C     (PZ S/PZ PY)
               W(69)=X(3)*(RI(13)*XX32+RI(14)*ZZ32)
     1           +RI(15)*(         +Z(3)*XZ32)
C     (PZ S/PZ PZ)
               W(70)=X(3)*(RI(13)*XX33+RI(14)*ZZ33)
     1           +RI(15)*(         +Z(3)*XZ33)
C     (PZ PX/S S)
               W(71)=RI(3)*XX31+RI(4)*ZZ31
C     (PZ PX/PX S)
               W(72)=X(1)*(RI(8)*XX31+RI(9)*ZZ31)
     1           +RI(10)*(Y(1)*XY31+Z(1)*XZ31)
C      (PZ PX/PX PX)
               W(73) =
     1     (RI(16)*XX31+RI(17)*ZZ31)*XX11+RI(18)*XX31*YYZZ11
     2     +RI(19)*ZZ31*ZZ11
     3     +RI(20)*(XY31*XY11+XZ31*XZ11)
     4     +RI(21)*ZZ31*YY11
     5     +RI(22)*YZ31*YZ11
C     (PZ PX/PY S)
               W(74)=X(2)*(RI(8)*XX31+RI(9)*ZZ31)
     1           +RI(10)*(Y(2)*XY31+Z(2)*XZ31)
C      (PZ PX/PY PX)
               W(75) =
     1     (RI(16)*XX31+RI(17)*ZZ31)*XX21+RI(18)*XX31*YYZZ21
     2     +RI(19)*ZZ31*ZZ21
     3     +RI(20)*(XY31*XY21+XZ31*XZ21)
     4     +RI(21)*ZZ31*YY21
     5     +RI(22)*YZ31*YZ21
C      (PZ PX/PY PY)
               W(76) =
     1     (RI(16)*XX31+RI(17)*ZZ31)*XX22+RI(18)*XX31*YYZZ22
     2     +RI(19)*ZZ31*ZZ22
     3     +RI(20)*(XY31*XY22+XZ31*XZ22)
     4     +RI(21)*ZZ31*YY22
     5     +RI(22)*YZ31*YZ22
C     (PZ PX/PZ S)
               W(77)=X(3)*(RI(8)*XX31+RI(9)*ZZ31)
     1           +RI(10)*(         +Z(3)*XZ31)
C     (PZ PX/PZ PX)
               W(78) =
     1      (RI(16)*XX31+RI(17)*ZZ31)*XX31
     2     +(RI(18)*XX31+RI(19)*ZZ31)*ZZ31
     3     +RI(20)*(XY31*XY31+XZ31*XZ31)
     4     +RI(22)*YZ31*YZ31
C      (PZ PX/PZ PY)
               W(79) =
     1      (RI(16)*XX31+RI(17)*ZZ31)*XX32
     2     +(RI(18)*XX31+RI(19)*ZZ31)*ZZ32
     3     +RI(20)*(XY31*XY32+XZ31*XZ32)
     4     +RI(22)*YZ31*YZ32
C      (PZ PX/PZ PZ)
               W(80) =
     1      (RI(16)*XX31+RI(17)*ZZ31)*XX33
     2     +(RI(18)*XX31+RI(19)*ZZ31)*ZZ33
     3     +RI(20)*XZ31*XZ33
C     (PZ PY/S S)
               W(81)=RI(3)*XX32+RI(4)*ZZ32
C     (PZ PY/PX S)
               W(82)=X(1)*(RI(8)*XX32+RI(9)*ZZ32)
     1           +RI(10)*(Y(1)*XY32+Z(1)*XZ32)
C      (PZ PY/PX PX)
               W(83) =
     1     (RI(16)*XX32+RI(17)*ZZ32)*XX11+RI(18)*XX32*YYZZ11
     2     +RI(19)*ZZ32*ZZ11
     3     +RI(20)*(XY32*XY11+XZ32*XZ11)
     4     +RI(21)*ZZ32*YY11
     5     +RI(22)*YZ32*YZ11
C     (PZ PY/PY S)
               W(84)=X(2)*(RI(8)*XX32+RI(9)*ZZ32)
     1           +RI(10)*(Y(2)*XY32+Z(2)*XZ32)
C      (PZ PY/PY PX)
               W(85) =
     1     (RI(16)*XX32+RI(17)*ZZ32)*XX21+RI(18)*XX32*YYZZ21
     2     +RI(19)*ZZ32*ZZ21
     3     +RI(20)*(XY32*XY21+XZ32*XZ21)
     4     +RI(21)*ZZ32*YY21
     5     +RI(22)*YZ32*YZ21
C      (PZ PY/PY PY)
               W(86) =
     1     (RI(16)*XX32+RI(17)*ZZ32)*XX22+RI(18)*XX32*YYZZ22
     2     +RI(19)*ZZ32*ZZ22
     3     +RI(20)*(XY32*XY22+XZ32*XZ22)
     4     +RI(21)*ZZ32*YY22
     5     +RI(22)*YZ32*YZ22
C     (PZ PY/PZ S)
               W(87)=X(3)*(RI(8)*XX32+RI(9)*ZZ32)
     1           +RI(10)*(         +Z(3)*XZ32)
C      (PZ PY/PZ PX)
               W(88) =
     1      (RI(16)*XX32+RI(17)*ZZ32)*XX31
     2     +(RI(18)*XX32+RI(19)*ZZ32)*ZZ31
     3     +RI(20)*(XY32*XY31+XZ32*XZ31)
     4     +RI(22)*YZ32*YZ31
C      (PZ PY/PZ PY)
               W(89) =
     1      (RI(16)*XX32+RI(17)*ZZ32)*XX32
     2     +(RI(18)*XX32+RI(19)*ZZ32)*ZZ32
     3     +RI(20)*(XY32*XY32+XZ32*XZ32)
     4     +RI(22)*YZ32*YZ32
C       (PZ PY/PZ PZ)
               W(90) =
     1      (RI(16)*XX32+RI(17)*ZZ32)*XX33
     2     +(RI(18)*XX32+RI(19)*ZZ32)*ZZ33
     3     +RI(20)*XZ32*XZ33
C     (PZ PZ/S S)
               W(91)=RI(3)*XX33+RI(4)*ZZ33
C     (PZ PZ/PX S)
               W(92)=X(1)*(RI(8)*XX33+RI(9)*ZZ33)
     1           +RI(10)*(          Z(1)*XZ33)
C       (PZ PZ/PX PX)
               W(93) =
     1     (RI(16)*XX33+RI(17)*ZZ33)*XX11+RI(18)*XX33*YYZZ11
     2     +RI(19)*ZZ33*ZZ11
     3     +RI(20)*XZ33*XZ11
     4     +RI(21)*ZZ33*YY11
C     (PZ PZ/PY S)
               W(94)=X(2)*(RI(8)*XX33+RI(9)*ZZ33)
     1           +RI(10)*(         +Z(2)*XZ33)
C       (PZ PZ/PY PX)
               W(95) =
     1     (RI(16)*XX33+RI(17)*ZZ33)*XX21+RI(18)*XX33*YYZZ21
     2     +RI(19)*ZZ33*ZZ21
     3     +RI(20)*XZ33*XZ21
     4     +RI(21)*ZZ33*YY21
C       (PZ PZ/PY PY)
               W(96) =
     1     (RI(16)*XX33+RI(17)*ZZ33)*XX22+RI(18)*XX33*YYZZ22
     2     +RI(19)*ZZ33*ZZ22
     3     +RI(20)*XZ33*XZ22
     4     +RI(21)*ZZ33*YY22
C     (PZ PZ/PZ S)
               W(97)=X(3)*(RI(8)*XX33+RI(9)*ZZ33)
     1           +RI(10)*(         +Z(3)*XZ33)
C       (PZ PZ/PZ PX)
               W(98) =
     1      (RI(16)*XX33+RI(17)*ZZ33)*XX31
     2     +(RI(18)*XX33+RI(19)*ZZ33)*ZZ31
     3     +RI(20)*XZ33*XZ31
C       (PZ PZ/PZ PY)
               W(99) =
     1      (RI(16)*XX33+RI(17)*ZZ33)*XX32
     2     +(RI(18)*XX33+RI(19)*ZZ33)*ZZ32
     3     +RI(20)*XZ33*XZ32
C       (PZ PZ/PZ PZ)
               W(100) =
     1      (RI(16)*XX33+RI(17)*ZZ33)*XX33
     2     +(RI(18)*XX33+RI(19)*ZZ33)*ZZ33
     3     +RI(20)*XZ33*XZ33
               KI = 100
            ELSE
C     (PX S/S S)
               W(2)=RI(2)*X(1)
C     (PX PX/S S)
               W(3)=RI(3)*XX11+RI(4)*YYZZ11
C     (PY S/S S)
               W(4)=RI(2)*X(2)
C     (PY PX/S S)
               W(5)=RI(3)*XX21+RI(4)*YYZZ21
C     (PY PY/S S)
               W(6)=RI(3)*XX22+RI(4)*YYZZ22
C     (PZ S/SS)
               W(7)=RI(2)*X(3)
C     (PZ PX/S S)
               W(8)=RI(3)*XX31+RI(4)*ZZ31
C     (PZ PY/S S)
               W(9)=RI(3)*XX32+RI(4)*ZZ32
C     (PZ PZ/S S)
               W(10)=RI(3)*XX33+RI(4)*ZZ33
               KI = 10
            END IF
         END IF
C
C *** NOW ROTATE THE NUCLEAR ATTRACTION INTEGRALS.
C *** THE STORAGE OF THE NUCLEAR ATTRACTION INTEGRALS  CORE(KL/IJ) IS
C     (SS/)=1,   (SO/)=2,   (OO/)=3,   (PP/)=4
C
         E1B(1)=-CSS1
         IF(NATORB(NI).EQ.4) THEN
            E1B(2) = -CSP1 *X(1)
            E1B(3) = -CPPS1*XX11-CPPP1*YYZZ11
            E1B(4) = -CSP1 *X(2)
            E1B(5) = -CPPS1*XX21-CPPP1*YYZZ21
            E1B(6) = -CPPS1*XX22-CPPP1*YYZZ22
            E1B(7) = -CSP1 *X(3)
            E1B(8) = -CPPS1*XX31-CPPP1*ZZ31
            E1B(9) = -CPPS1*XX32-CPPP1*ZZ32
            E1B(10)= -CPPS1*XX33-CPPP1*ZZ33
         END IF
         E2A(1)=-CSS2
         IF(NATORB(NJ).EQ.4) THEN
            E2A(2) = -CSP2 *X(1)
            E2A(3) = -CPPS2*XX11-CPPP2*YYZZ11
            E2A(4) = -CSP2 *X(2)
            E2A(5) = -CPPS2*XX21-CPPP2*YYZZ21
            E2A(6) = -CPPS2*XX22-CPPP2*YYZZ22
            E2A(7) = -CSP2 *X(3)
            E2A(8) = -CPPS2*XX31-CPPP2*ZZ31
            E2A(9) = -CPPS2*XX32-CPPP2*ZZ32
            E2A(10)= -CPPS2*XX33-CPPP2*ZZ33
         END IF
         IF(ABS(TORE(NI)).GT.20.D0.AND.ABS(TORE(NJ)).GT.20.D0) THEN
C SPARKLE-SPARKLE INTERACTION
            ENUC=0.D0
            RETURN
         ELSEIF (RIJ.LT.1.D0.AND.NATORB(NI)*NATORB(NJ).EQ.0) THEN
            ENUC=0.D0
            RETURN
         ENDIF
         SCALE = EXP(-ALP(NI)*RIJ)+EXP(-ALP(NJ)*RIJ)
C
         IF (NI.EQ.24.AND.NJ.EQ.24) THEN
            SCALE = EXP(-ALPTM(NI)*RIJ)+EXP(-ALPTM(NJ)*RIJ)
         ENDIF
C
         NT=NI+NJ
         IF(NT.EQ.8.OR.NT.EQ.9) THEN
            IF(NI.EQ.7.OR.NI.EQ.8) SCALE=SCALE+(RIJ-1.D0)*EXP(-ALP(NI)*R
     1IJ)
            IF(NJ.EQ.7.OR.NJ.EQ.8) SCALE=SCALE+(RIJ-1.D0)*EXP(-ALP(NJ)*R
     1IJ)
         ENDIF
         ENUC = TORE(NI)*TORE(NJ)*GAM
         SCALE=ABS(SCALE*ENUC)
         IF( AM1PM3 )THEN
         IF(AM1.AND.(NI.EQ.5.OR.NJ.EQ.5))THEN
C
C   LOAD IN AM1 BORON GAUSSIANS
C
          NK=NI+NJ-5
C   NK IS THE ATOMIC NUMBER OF THE NON-BORON ATOM
          NL=1
          IF(NK.EQ.1)NL=2
          IF(NK.EQ.6)NL=3
          IF(NK.EQ.9.OR.NK.EQ.17.OR.NK.EQ.35.OR.NK.EQ.53)NL=4
          DO 25 I=1,3
          FN1(5,I)=BORON1(I,NL)
          FN2(5,I)=BORON2(I,NL)
  25      FN3(5,I)=BORON3(I,NL)
          ENDIF
            DO 30 IG=1,10
               IF(ABS(FN1(NI,IG)).GT.0.D0) THEN
                  AX = FN2(NI,IG)*(RIJ-FN3(NI,IG))**2
                  IF(AX .LE. 25.D0) THEN
                     SCALE=SCALE +TORE(NI)*TORE(NJ)/RIJ*FN1(NI,IG)*EXP(-
     1AX)
                  ENDIF
               ENDIF
               IF(ABS(FN1(NJ,IG)).GT.0.D0) THEN
                  AX = FN2(NJ,IG)*(RIJ-FN3(NJ,IG))**2
                  IF(AX .LE. 25.D0) THEN
                     SCALE=SCALE +TORE(NI)*TORE(NJ)/RIJ*FN1(NJ,IG)*EXP(-
     1AX)
                  ENDIF
               ENDIF
   30       CONTINUE
         ENDIF
         ENUC=ENUC+SCALE
C
         IF(NATORB(NI)*NATORB(NJ).EQ.0)KI=0
         KR=KR+KI
C
C
      ENDIF
      RETURN
      END
      SUBROUTINE RSP(A,N,MATZ,W,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION A(*),  W(N), Z(N,N)
*******************************************************************
*
*   EISPACK DIAGONALIZATION ROUTINES: TO FIND THE EIGENVALUES AND
*           EIGENVECTORS (IF DESIRED) OF A REAL SYMMETRIC PACKED MATRIX.
* ON INPUT-      N  IS THE ORDER OF THE MATRIX  A,
*                A  CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
*                   PACKED MATRIX STORED ROW-WISE,
*             MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF ONLY
*                   EIGENVALUES ARE DESIRED,  OTHERWISE IT IS SET TO
*                   ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND
*                   EIGENVECTORS.
* ON OUTPUT-     W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER,
*                Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO,
*
*******************************************************************
* THIS SUBROUTINE WAS CHOSEN AS BEING THE MOST RELIABLE. (JJPS)
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C FV and FV2 are scratch arrays which shouldo be dimensioned to the 
C to the order of the biggest diagonalizable matrix
C This is accomplished by the parameter MAXDMO ( Max Diag. Mat. order)   IR0494
      DIMENSION FV1(MAXDMO),FV2(MAXDMO)
        SAVE
C   Test that MAXDMO is big enough                                       IR0494
      MX=MAX(MAXORB,NMECI*NMECI)
      IF(MAXDMO.LT.MX) then
        WRITE(6,'(3x,"*** MAXDMO should be at least",I5," ***")') MX
        WRITE(6,'("Please, change MAXDMO in SIZES.i and recompile")')
        STOP
      ENDIF
      NV=(MAXDMO*(MAXDMO+1))/2                                           IR0494
      NM=N
C
20    call tred3(n,nv,a,w,fv1,fv2)
      IF (MATZ .NE. 0) GO TO 30
C     ********** FIND EIGENVALUES ONLY **********
      call tqlrat(n,w,fv2,ierr)
      GO TO 60
C     ********** FIND BOTH EIGENVALUES AND EIGENVECTORS **********
   30 DO 50    I = 1, N
C
         DO 40    J = 1, N
            Z(J,I)=0.0D0
   40    CONTINUE
C
         Z(I,I)=1.0D0
   50 CONTINUE
C
      call tql2(nm,n,w,fv1,z,ierr)
      IF (IERR .NE. 0) GO TO 70
      call trbak3(nm,n,nv,a,n,z,ierr)                                   IR0494
   60 return
C     ********** LAST CARD OF RSP **********
   70 write(6,'( "Error",I3," in RSP" )') ierr                          IR0494
      RETURN
      END
      SUBROUTINE SCHMIB(U,N,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
C
C     SAME AS SCHMIDT BUT WORKS FROM RIGHT TO LEFT.
C
      DIMENSION U(NDIM,NDIM)
         SAVE                                                           GL0892
      DATA ZERO,SMALL,ONE/0.0,0.01,1.0/
      N1=N+1
      II=0
      DO 110 K=1,N
         K1=K-1
C
C     NORMALIZE KTH COLUMN VECTOR
C
         DOT = ZERO
         DO 10 I=1,N
   10    DOT=DOT+U(I,N1-K)*U(I,N1-K)
         IF(DOT.EQ.ZERO) GO TO 100
         SCALE=ONE/SQRT(DOT)
         DO 20 I=1,N
   20    U(I,N1-K)=SCALE*U(I,N1-K)
   30    IF(K1.EQ.0) GO TO 110
         NPASS=0
C
C     PROJECT OUT K-1 PREVIOUS ORTHONORMAL VECTORS FROM KTH VECTOR
C
   40    NPASS=NPASS+1
         DO 70 J=1,K1
            DOT=ZERO
            DO 50 I=1,N
   50       DOT=DOT+U(I,N1-J)*U(I,N1-K)
            DO 60 I=1,N
   60       U(I,N1-K)=U(I,N1-K)-DOT*U(I,N1-J)
   70    CONTINUE
C
C     SECOND NORMALIZATION (AFTER PROJECTION)
C     IF KTH VECTOR IS SMALL BUT NOT ZERO THEN NORMALIZE
C     AND PROJECT AGAIN TO CONTROL ROUND-OFF ERRORS.
C
         DOT=ZERO
         DO 80 I=1,N
   80    DOT=DOT+U(I,N1-K)*U(I,N1-K)
         IF(DOT.EQ.ZERO) GO TO 100
         IF(DOT.LT.SMALL.AND.NPASS.GT.2) GO TO 100
         SCALE=ONE/SQRT(DOT)
         DO 90 I=1,N
   90    U(I,N1-K)=SCALE*U(I,N1-K)
         IF(DOT.LT.SMALL) GO TO 40
         GO TO 110
C
C     REPLACE LINEARLY DEPENDENT KTH VECTOR BY A UNIT VECTOR.
C
  100    II=II+1
C     IF(II.GT.N) STOP
         U(II,N1-K)=ONE
         GO TO 30
  110 CONTINUE
      RETURN
      END
      SUBROUTINE SCHMIT(U,N,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION U(NDIM,NDIM)
         SAVE                                                           GL0892
      DATA ZERO,SMALL,ONE/0.0,0.01,1.0/
      II=0
      DO 110 K=1,N
         K1=K-1
C
C     NORMALIZE KTH COLUMN VECTOR
C
         DOT = ZERO
         DO 10 I=1,N
   10    DOT=DOT+U(I,K)*U(I,K)
         IF(DOT.EQ.ZERO) GO TO 100
         SCALE=ONE/SQRT(DOT)
         DO 20 I=1,N
   20    U(I,K)=SCALE*U(I,K)
   30    IF(K1.EQ.0) GO TO 110
         NPASS=0
C
C     PROJECT OUT K-1 PREVIOUS ORTHONORMAL VECTORS FROM KTH VECTOR
C
   40    NPASS=NPASS+1
         DO 70 J=1,K1
            DOT=ZERO
            DO 50 I=1,N
   50       DOT=DOT+U(I,J)*U(I,K)
            DO 60 I=1,N
   60       U(I,K)=U(I,K)-DOT*U(I,J)
   70    CONTINUE
C
C     SECOND NORMALIZATION (AFTER PROJECTION)
C     IF KTH VECTOR IS SMALL BUT NOT ZERO THEN NORMALIZE
C     AND PROJECT AGAIN TO CONTROL ROUND-OFF ERRORS.
C
         DOT=ZERO
         DO 80 I=1,N
   80    DOT=DOT+U(I,K)*U(I,K)
         IF(DOT.EQ.ZERO) GO TO 100
         IF(DOT.LT.SMALL.AND.NPASS.GT.2) GO TO 100
         SCALE=ONE/SQRT(DOT)
         DO 90 I=1,N
   90    U(I,K)=SCALE*U(I,K)
         IF(DOT.LT.SMALL) GO TO 40
         GO TO 110
C
C     REPLACE LINEARLY DEPENDENT KTH VECTOR BY A UNIT VECTOR.
C
  100    II=II+1
C     IF(II.GT.N) STOP
         U(II,K)=ONE
         GO TO 30
  110 CONTINUE
      RETURN
      END
      SUBROUTINE SET (S1,S2,NA,NB,RAB,NBOND,II)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /SETC/ A(7),B(7),SA,SB,FACTOR,ISP,IPS
      DIMENSION TEMP(20)
      DIMENSION FACT(17)
         SAVE                                                           GL0892
      DATA FACT/1.D0,2.D0,6.D0,24.D0,120.D0,720.D0,5040.D0,40320.D0,
     1362880.D0,3628800.D0,39916800.D0,479001600.D0,6227020800.D0,
     28.71782912D10,1.307674368D12,2.092278989D13,3.556874281D14/
C***********************************************************************
C
C     SET IS PART OF THE OVERLAP CALCULATION, CALLED BY OVERLP.
C
C***********************************************************************
      IF (NA.GT.NB) then
        ISP=2
        IPS=1
        SA=S2
        SB=S1
      else
        ISP=1
        IPS=2
        SA=S1
        SB=S2
      endif
      J=II+2
      IF (II.GT.3) J=J-1
      ALPHA=0.5D00*RAB*(SA+SB)
      BETA=0.5D00*RAB*(SB-SA)
      JCALL=J-1
      c = exp(-alpha)
      a(1) = c/alpha
      do 500 i=1,jcall
      a(i+1) = (a(i)*i+c)/alpha
500   continue
c------------------
c---- call BINTGS(beta,jcall)
c------------------
      ABSX = ABS(beta)
      IF (ABSX.GT.3.D00) GO TO 40
      IF (ABSX.LE.2.D00) GO TO 10
      IF (jcall.LE.10) GO TO 40
      LAST=15
      GO TO 60
   10 IF (ABSX.LE.1.D00) GO TO 20
      IF (jcall.LE.7) GO TO 40
      LAST=12
      GO TO 60
   20 IF (ABSX.LE.0.5D00) GO TO 30
      IF (jcall.le.5) GO TO 40
      LAST=7
      GO TO 60
   30 IF (ABSX.LE.1.D-6) GOTO 90
      LAST=6
      GO TO 60
   40 EXPX=EXP(beta)
      EXPMX=1.D00/EXPX
      xinv = 1.0/beta
      b(1) = (expx - expmx)*xinv
      zk = 1.0
      DO 50 I=1,jcall
      zk = -1.0*zk
      b(i+1) = (i*b(i) + zk*expx - expmx)*xinv
50    continue
      return
   60 DO 80 I=0,jcall
          zk = -beta
          y = 2.0*mod(i+1,2)/(i+1)
          do 70 m=1,last
          y = y + 2.0*zk*mod(m+i+1,2)/(fact(m)*(m+i+1))
          zk = -beta*zk
70        continue
   80 B(I+1)=Y
      return
   90 DO 100 I=0,jcall
  100 B(I+1)=(2*MOD(I+1,2))/(I+1.D0)
c--------------
      RETURN
C
      END
      SUBROUTINE SETUPG
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      COMMON /EXPONT/ ZS(107),ZP(107),ZD(107)
      COMMON /NATYPE/ NZTYPE(107),MTYPE(30),LTYPE
      COMMON /TEMP/  CC(60,6),ZZ(60,6)
      DIMENSION ALLC(6,5,2), ALLZ(6,5,2)
         SAVE                                                           GL0892
C     SET-UP THE STEWART'S STO-6G EXPANSIONS
C                                            1S
      ALLZ(1,1,1) =2.310303149D01
      ALLZ(2,1,1) =4.235915534D00
      ALLZ(3,1,1) =1.185056519D00
      ALLZ(4,1,1) =4.070988982D-01
      ALLZ(5,1,1) =1.580884151D-01
      ALLZ(6,1,1) =6.510953954D-02
C
      ALLC(1,1,1) =9.163596280D-03
      ALLC(2,1,1) =4.936149294D-02
      ALLC(3,1,1) =1.685383049D-01
      ALLC(4,1,1) =3.705627997D-01
      ALLC(5,1,1) =4.164915298D-01
      ALLC(6,1,1) =1.303340841D-01
C                                      2S
      ALLZ(1,2,1) =2.768496241D01
      ALLZ(2,2,1) =5.077140627D00
      ALLZ(3,2,1) =1.426786050D00
      ALLZ(4,2,1) =2.040335729D-01
      ALLZ(5,2,1) =9.260298399D-02
      ALLZ(6,2,1) =4.416183978D-02
C
      ALLC(1,2,1) =-4.151277819D-03
      ALLC(2,2,1) =-2.067024148D-02
      ALLC(3,2,1) =-5.150303337D-02
      ALLC(4,2,1) =3.346271174D-01
      ALLC(5,2,1) =5.621061301D-01
      ALLC(6,2,1) =1.712994697D-01
C                                     2P
      ALLZ(1,2,2) =5.868285913D00
      ALLZ(2,2,2) =1.530329631D00
      ALLZ(3,2,2) =5.475665231D-01
      ALLZ(4,2,2) =2.288932733D-01
      ALLZ(5,2,2) =1.046655969D-01
      ALLZ(6,2,2) =4.948220127D-02
C
      ALLC(1,2,2) =7.924233646D-03
      ALLC(2,2,2) =5.144104825D-02
      ALLC(3,2,2) =1.898400060D-01
      ALLC(4,2,2) =4.049863191D-01
      ALLC(5,2,2) =4.012362861D-01
      ALLC(6,2,2) =1.051855189D-01
C                                      3S
      ALLZ(1,3,1) =3.273031938D00
      ALLZ(2,3,1) =9.200611311D-01
      ALLZ(3,3,1) =3.593349765D-01
      ALLZ(4,3,1) =8.636686991D-02
      ALLZ(5,3,1) =4.797373812D-02
      ALLZ(6,3,1) =2.724741144D-02
      ALLC(1,3,1) =-6.775596947D-03
      ALLC(2,3,1) =-5.639325779D-02
      ALLC(3,3,1) =-1.587856086D-01
      ALLC(4,3,1) =5.534527651D-01
      ALLC(5,3,1) =5.015351020D-01
      ALLC(6,3,1) =7.223633674D-02
C                                     3P
      ALLZ(1,3,2) =5.077973607D00
      ALLZ(2,3,2) =1.340786940D00
      ALLZ(3,3,2) =2.248434849D-01
      ALLZ(4,3,2) =1.131741848D-01
      ALLZ(5,3,2) =6.076408893D-02
      ALLZ(6,3,2) =3.315424265D-02
      ALLC(1,3,2) =-3.329929840D-03
      ALLC(2,3,2) =-1.419488340D-02
      ALLC(3,3,2) =1.639395770D-01
      ALLC(4,3,2) =4.485358256D-01
      ALLC(5,3,2) =3.908813050D-01
      ALLC(6,3,2) =7.411456232D-02
C                                     4S
      ALLZ(1,4,1) = 1.365346 D+00
      ALLZ(2,4,1) = 4.393213 D-01
      ALLZ(3,4,1) = 1.877069 D-01
      ALLZ(4,4,1) = 9.360270 D-02
      ALLZ(5,4,1) = 5.052263 D-02
      ALLZ(6,4,1) = 2.809354 D-02
      ALLC(1,4,1) = 3.775056 D-03
      ALLC(2,4,1) =-5.585965 D-02
      ALLC(3,4,1) =-3.192946 D-01
      ALLC(4,4,1) =-2.764780 D-02
      ALLC(5,4,1) = 9.049199 D-01
      ALLC(6,4,1) = 3.406258 D-01
C                                   4P
      ALLC(1,4,2) =-7.052075 D-03
      ALLC(2,4,2) =-5.259505 D-02
      ALLC(3,4,2) =-3.773450 D-02
      ALLC(4,4,2) = 3.874773 D-01
      ALLC(5,4,2) = 5.791672 D-01
      ALLC(6,4,2) = 1.221817 D-01
      ALLZ(1,4,2) = 1.365346 D+00
      ALLZ(2,4,2) = 4.393213 D-01
      ALLZ(3,4,2) = 1.877069 D-01
      ALLZ(4,4,2) = 9.360270 D-02
      ALLZ(5,4,2) = 5.052263 D-02
      ALLZ(6,4,2) = 2.809354 D-02
C                                     5S
      ALLZ(1,5,1) = 7.701420258D-01
      ALLZ(2,5,1) = 2.756268915D-01
      ALLZ(3,5,1) = 1.301847480D-01
      ALLZ(4,5,1) = 6.953441940D-02
      ALLZ(5,5,1) = 4.002545502D-02
      ALLZ(6,5,1) = 2.348388309D-02
      ALLC(1,5,1) = 1.267447151D-02
      ALLC(2,5,1) = 3.266734789D-03
      ALLC(3,5,1) =-4.307553999D-01
      ALLC(4,5,1) =-3.231998963D-01
      ALLC(5,5,1) = 1.104322879D+00
      ALLC(6,5,1) = 4.368498703D-01
C                                      5P
      ALLZ(1,5,2) = 7.701420258D-01
      ALLZ(2,5,2) = 2.756268915D-01
      ALLZ(3,5,2) = 1.301847480D-01
      ALLZ(4,5,2) = 6.953441940D-02
      ALLZ(5,5,2) = 4.002545502D-02
      ALLZ(6,5,2) = 2.348388309D-02
      ALLC(1,5,2) =-1.105673292D-03
      ALLC(2,5,2) =-6.243132446D-02
      ALLC(3,5,2) =-1.628476766D-01
      ALLC(4,5,2) = 3.210328714D-01
      ALLC(5,5,2) = 6.964579592D-01
      ALLC(6,5,2) = 1.493146125D-01
      DO 30 I=1,10
         IF(MTYPE(I).EQ.0)GOTO 30
         NI=MTYPE(I)
         XI=ZS(NI)
         IA=I*4-3
         IB=IA+3
         IF(NI.LT.2) THEN
            NQN=1
         ELSEIF(NI.LT.10)THEN
            NQN=2
         ELSEIF(NI.LT.18)THEN
            NQN=3
         ELSEIF(NI.LT.36)THEN
            NQN=4
         ELSEIF(NI.LT.54)THEN
            NQN=5
         ELSE
            WRITE(6,*)' NO GAUSSIANS AVAILABLE'
            STOP
         ENDIF
         DO 20 K=IA,IB
            L=1
            IF(K.GT.IA) L=2
            IF(K.GT.IA) XI=ZP(NI)
            DO 10 J=1,6
               CC(K,J)=ALLC(J,NQN,L)
   10       ZZ(K,J)=ALLZ(J,NQN,L)*XI**2
   20    CONTINUE
   30 CONTINUE
      RETURN
      END
      FUNCTION SPCG(C1,C2,C3,C4,W,WJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION C1(*),C2(*),C3(*),C4(*),W(*), WJ(*)
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /TWOELE/ GSS(107),GSP(107),GPP(107),GP2(107),HSP(107)
     1                ,GSD(107),GPD(107),GDD(107)
      COMMON /EULER / TVEC(3,3), ID
C********************************************************************
C
C     SPCG CALCULATES THE REPULSION BETWEEN ELECTRON 1 IN MOLECULAR
C     ORBITALS C1 AND C2 AND ELECTRON 2 IN M.O.S C3 AND C4 FOR THE
C     VALENCE SP SHELL AT AN MNDO OR MINDO/3 LEVEL.
C
C                            USAGE
C      XJ=SPCG(C(1,I),C(1,J),C(1,K),C(1,L))
C  OR, XJ=<I(1),J(1)/K(2),L(2)>
C
C    ON INPUT C1    THE FIRST COLUMN MOLECULAR ORBITAL OF ELECTRON ONE.
C             C2        SECOND
C             C3        FIRST                                      TWO.
C             C4        SECOND
C
C   ON OUTPUT SPCG   =   <C1(1)*C2(1)/C3(2)*C4(2)>
C*********************************************************************
      COMMON /KEYWRD/ KEYWRD
      CHARACTER*80 KEYWRD
      LOGICAL LID
         SAVE                                                           GL0892
      DATA ITYPE /1/
   10 CONTINUE
      GOTO (20,70,30) ITYPE
   20 CONTINUE
      LID=(ID.EQ.0)
      IF(INDEX(KEYWRD,'MINDO') .NE. 0) THEN
         ITYPE=2
      ELSE
         ITYPE=3
      ENDIF
      GOTO 10
*                           ******************
*                           *      MNDO      *
*                           *     OPTION     *
*                           ******************
   30 SPCG=0.0D00
      KK=0
      DO 60 II=1,NUMAT
         IA=NFIRST(II)
         IB=NLAST(II)
         IMINUS=II-1
         DO 50 JJ=1,IMINUS
            JA=NFIRST(JJ)
            JB=NLAST(JJ)
            DO 40 I=IA,IB
               DO 40 J=IA,I
                  DO 40 K=JA,JB
                     DO 40 L=JA,K
                        KK=KK+1
                        IF(LID) THEN
                           WINT=W(KK)
                        ELSE
                           WINT=WJ(KK)
                        ENDIF
                        SPCG=SPCG+WINT*(C1(I)*C2(J)*C3(K)*C4(L)
     1 + C1(K)*C2(L)*C3(I)*C4(J))
                        IF(I.NE.J)
     1SPCG=SPCG+WINT*(C1(J)*C2(I)*C3(K)*C4(L)
     2 + C1(K)*C2(L)*C3(J)*C4(I))
                        IF(K.NE.L) SPCG=SPCG+WINT*(C1(I)*C2(J)*C3(L)*C4(
     1K) + C1(L)*C2(K)*C3(I)*C4(J))
                        IF((I.NE.J).AND.(K.NE.L))SPCG=SPCG+WINT*(C1(J)*C
     12(I)*C3(L)*C4(K) +C1(L)*C2(K)*C3(J)*C4(I))
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
      GOTO 110
*                           ******************
*                           *     MINDO/3    *
*                           *     OPTION     *
*                           ******************
   70 CONTINUE
      SPCG=0.D0
      KR=0
      DO 100 II=1,NUMAT
         IA=NFIRST(II)
         IB=NLAST(II)
         IM1=II-1
         DO 90 JJ=1,IM1
            KR=KR+1
            IF(LID) THEN
               ELREP=W(KR)
            ELSE
               ELREP=WJ(KR)
            ENDIF
            JA=NFIRST(JJ)
            JB=NLAST(JJ)
            DO 80 I=IA,IB
               DO 80 K=JA,JB
   80       SPCG=SPCG+ELREP*(C1(I)*C2(I)*C3(K)*C4(K)+
     1                                 C1(K)*C2(K)*C3(I)*C4(I))
   90    CONTINUE
  100 CONTINUE
  110 CONTINUE
      ATEMP=SPCG
      IS1=0
      DO 150 I1=1,NUMAT
         IS1=IS1+1
         IZN=NAT(I1)
C
C      (SS/SS)
C
         SPCG=SPCG+C1(IS1)*C2(IS1)*C3(IS1)*C4(IS1)*GSS(IZN)
         IF(IZN.LT.3) GO TO 150
         IS=IS1
         IS1=IS1+1
         IX=IS1
         IS1=IS1+1
         IY=IS1
         IS1=IS1+1
         IZ=IS1
         SPCG=SPCG+GPP(IZN)*
     1                       (
     2           C1(IX)*C2(IX)*C3(IX)*C4(IX)+
     3           C1(IY)*C2(IY)*C3(IY)*C4(IY)+
     4           C1(IZ)*C2(IZ)*C3(IZ)*C4(IZ)
     5                       )
         SPCG=SPCG+GSP(IZN)*
     1                       (
     2           C1(IS)*C2(IS)*C3(IX)*C4(IX)+
     3           C1(IS)*C2(IS)*C3(IY)*C4(IY)+
     4           C1(IS)*C2(IS)*C3(IZ)*C4(IZ)+
     5           C1(IX)*C2(IX)*C3(IS)*C4(IS)+
     6           C1(IY)*C2(IY)*C3(IS)*C4(IS)+
     7           C1(IZ)*C2(IZ)*C3(IS)*C4(IS)
     8                       )
         SPCG=SPCG+GP2(IZN)*
     1                       (
     2           C1(IX)*C2(IX)*C3(IY)*C4(IY)+
     3           C1(IX)*C2(IX)*C3(IZ)*C4(IZ)+
     4           C1(IY)*C2(IY)*C3(IZ)*C4(IZ)+
     5           C1(IY)*C2(IY)*C3(IX)*C4(IX)+
     6           C1(IZ)*C2(IZ)*C3(IX)*C4(IX)+
     7           C1(IZ)*C2(IZ)*C3(IY)*C4(IY)
     8                       )
         TEMP1=HSP(IZN)
         DO 120 J1=IX,IZ
            SPCG=SPCG+TEMP1*
     1                       (
     2           C1(IS)*C2(J1)*C3(J1)*C4(IS)+
     3           C1(IS)*C2(J1)*C3(IS)*C4(J1)+
     4           C1(J1)*C2(IS)*C3(IS)*C4(J1)+
     5           C1(J1)*C2(IS)*C3(J1)*C4(IS)
     6                       )
  120    CONTINUE
         TEMP1=0.5D0*(GPP(IZN)-GP2(IZN))
         DO 140 J1=IX,IZ
            DO 130 K1=IX,IZ
               IF(J1.EQ.K1) GO TO 130
               SPCG=SPCG+TEMP1*
     1                       (
     2           C1(J1)*C2(K1)*C3(J1)*C4(K1)+
     3           C1(J1)*C2(K1)*C3(K1)*C4(J1)
     4                       )
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
      SUBROUTINE SPLINE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      LOGICAL SKIP1,SKIP2
C
C     FIT F(X) BY A CUBIC SPLINE GIVEN VALUES OF THE FUNCTION
C     AND ITS FIRST DERIVATIVE AT N PNTS.
C     SUBROUTINE RETURNS VALUES OF XMIN,FMIN, AND DFMIN
C     AND MAY REORDER THE DATA.
C     CALLING PROGRAM SUPPLIES ALL OTHER VALUES IN THE
C     COMMON BLOCK.
C     XLOW AND XHIGH SET LIMITS ON THE INTERVAL WITHIN WHICH
C     TO SEARCH.  SUBROUTINE MAY FURTHER REDUCE THIS INTERVAL.
C
C common FIT splitted for portability (Ivan Rossi 0394 &8) )
      COMMON /FITPTS/ N
     1       /FITDAT/ XLOW,XHIGH,XMIN,FMIN,DFMIN,X(12),F(12),DF(12)
         SAVE                                                           GL0892
      DATA CLOSE, BIG, HUGE, USTEP, DSTEP/1.0E-8,500.0,1.0E+10,1.0,2.0/
C
C     SUBROUTINE ASSUMES THAT THE FIRST N-1 DATA PNTS HAVE BEEN
C     PREVIOUSLY ORDERED,  X(I).LT.X(I+1) FOR I=1,2,...,N-2
C     NOW MOVE NTH POINT TO ITS PROPER PLACE.
C
      XMIN=X(N)
      FMIN=F(N)
      DFMIN=DF(N)
      N1=N-1
      K=N1
   10 IF(X(K).LT.XMIN) GO TO 20
      X(K+1)=X(K)
      F(K+1)=F(K)
      DF(K+1)=DF(K)
      K=K-1
      IF(K.GT.0) GO TO 10
   20 X(K+1)=XMIN
      F(K+1)=FMIN
      DF(K+1)=DFMIN
C
C     DEFINE THE INTERVAL WITHIN WHICH WE TRUST THE SPLINE FIT.
C     USTEP =  UP HILL STEP SIZE FACTOR
C     DSTEP = DOWN HILL STEP SIZE FACTOR
C
      IF(DF(1).GT.0.0) STEP=DSTEP
      IF(DF(1).LE.0.0) STEP=USTEP
      XSTART=X(1)-STEP*(X(2)-X(1))
      XSTART=MAX(XSTART,XLOW)
      IF(DF(N).GT.0.0) STEP=USTEP
      IF(DF(N).LE.0.0) STEP=DSTEP
      XSTOP=X(N)+STEP*(X(N)-X(N1))
      XSTOP=MIN(XSTOP,XHIGH)
C
C     SEARCH FOR MINIMUM
C
      DO 110 K=1,N1
         SKIP1=K.NE.1
         SKIP2=K.NE.N1
         IF(F(K).GE.FMIN) GO TO 30
         XMIN=X(K)
         FMIN=F(K)
         DFMIN=DF(K)
   30    DX=X(K+1)-X(K)
C
C     SKIP INTERVAL IF PNTS ARE TOO CLOSE TOGETHER
C
         IF(DX.LE.CLOSE) GO TO 110
         X1=0.0
         IF(K.EQ.1) X1=XSTART-X(1)
         X2=DX
         IF(K.EQ.N1) X2=XSTOP-X(N1)
C
C     (A,B,C)=COEF OF (CUBIC,QUADRATIC,LINEAR) TERMS
C
         DUM=(F(K+1)-F(K))/DX
         A=(DF(K)+DF(K+1)-DUM-DUM)/(DX*DX)
         B=(DUM+DUM+DUM-DF(K)-DF(K)-DF(K+1))/DX
         C=DF(K)
C
C     XK = X-X(K) AT THE MINIMUM WITHIN THE KTH SUBINTERVAL
C     TEST FOR PATHOLOGICAL CASES.
C
         BB=B*B
         AC3=(A+A+A)*C
         IF(BB.LT.AC3) GO TO 90
         IF( B.GT.0.0) GO TO 40
         IF(ABS(B).GT.HUGE*ABS(A)) GO TO 90
         GO TO 50
   40    IF(BB.GT.BIG*ABS(AC3)) GO TO 60
C
C     WELL BEHAVED CUBIC
C
   50    XK=(-B+SQRT(BB-AC3))/(A+A+A)
         GO TO 70
C
C     CUBIC IS DOMINATED BY QUADRATIC TERM
C
   60    R=AC3/BB
         XK=-(((0.039063*R+0.0625)*R+0.125)*R+0.5)*C/B
   70    IF(XK.LT.X1.OR.XK.GT.X2) GO TO 90
   80    FM=((A*XK+B)*XK+C)*XK+F(K)
         IF(FM.GT.FMIN) GO TO 90
         XMIN=XK+X(K)
         FMIN=FM
         DFMIN=((A+A+A)*XK+B+B)*XK+C
C
C     EXTRAPOLATE TO END OF INTERVAL IF K=1 AND/OR K=N1
C
   90    IF(SKIP1) GO TO 100
         SKIP1=.TRUE.
         XK=X1
         GO TO 80
  100    IF(SKIP2) GO TO 110
         SKIP2=.TRUE.
         XK=X2
         GO TO 80
  110 CONTINUE
      RETURN
      END
      SUBROUTINE SWAP(C,N,MDIM,NOCC,IFILL)
      IMPLICIT  DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION C(MDIM,MDIM)
C******************************************************************
C
C        SWAP ENSURES THAT A NAMED MOLECULAR ORBITAL IFILL IS FILLED
C ON INPUT
C          C = EIGENVECTORS IN A MDIM*MDIM MATRIX
C          N = NUMBER OF ORBITALS
C          NOCC = NUMBER OF OCCUPIED ORBITALS
C          IFILL = FILLED ORBITAL
C******************************************************************
      COMMON /SWAP0/ PSI(MAXORB), STDPSI(MAXORB)
         SAVE                                                           GL0892
      IF(IFILL.GT.0) GOTO 20
C
C     WE NOW DEFINE THE FILLED ORBITAL
C
      IFILL=-IFILL
      DO 10 I=1,N
         STDPSI(I)=C(I,IFILL)
   10 PSI(I)=C(I,IFILL)
      RETURN
   20 CONTINUE
C
C     FIRST FIND THE LOCATION OF IFILL
C
      SUM=0.D0
      DO 30 I=1,N
   30 SUM=SUM+PSI(I)*C(I,IFILL)
      IF(ABS(SUM).GT.0.7071D0) GOTO 90
C
C     IFILL HAS MOVED!
C
      SUMMAX=0.D0
      DO 50 IFILL=1,N
         SUM=0.D0
         DO 40 I=1,N
   40    SUM=SUM+STDPSI(I)*C(I,IFILL)
         SUM=ABS(SUM)
         IF(SUM.GT.SUMMAX)JFILL=IFILL
         IF(SUM.GT.SUMMAX)SUMMAX=SUM
         IF(SUM.GT.0.7071D0) GOTO 90
   50 CONTINUE
      DO 70 IFILL=1,N
         SUM=0.D0
         DO 60 I=1,N
   60    SUM=SUM+PSI(I)*C(I,IFILL)
         SUM=ABS(SUM)
         IF(SUM.GT.SUMMAX)JFILL=IFILL
         IF(SUM.GT.SUMMAX)SUMMAX=SUM
         IF(SUM.GT.0.7071D0) GOTO 90
   70 CONTINUE
      WRITE(6,80)SUMMAX,JFILL
   80 FORMAT(/,' CAUTION !!! SUM IN SWAP VERY SMALL, SUMMAX =',F10.5,
     1' JFILL=',I3)
      IFILL=JFILL
   90 CONTINUE
      IF(IFILL.LE.NOCC) RETURN
C
C    ITS EMPTY, SO SWAP IT WITH THE HIGHEST FILLED
C
      DO 100 I=1,N
         X=C(I,NOCC)
         C(I,NOCC)=C(I,IFILL)
         C(I,IFILL)=X
  100 CONTINUE
      RETURN
      END
      SUBROUTINE SYMTRY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /GEOSYM/ NDEP, LOCPAR(MAXPAR), IDEPFN(MAXPAR),
     1         LOCDEP(MAXPAR)
         SAVE                                                           GL0892
C**********************************************************************
C
C  SYMTRY COMPUTES THE BOND LENGTHS AND ANGLES THAT ARE FUNCTIONS OF
C         OTHER BOND LENGTHS AND ANGLES.
C
C ON INPUT GEO     = KNOWN INTERNAL COORDINATES
C          NDEP    = NUMBER OF DEPENDENCY FUNCTIONS.
C          IDEPFN  = ARRAY OF DEPENDENCY FUNCTIONS.
C          LOCDEP  = ARRAY OF LABELS OF DEPENDENT ATOMS.
C          LOCPAR  = ARRAY OF LABELS OF REFERENCE ATOMS.
C
C  ON OUTPUT THE ARRAY "GEO" IS FILLED
C***********************************************************************
C
C     NOW COMPUTE THE DEPENDENT PARAMETERS.
C
      DO 10 I=1,NDEP
         CALL HADDON (VALUE,LOCN,IDEPFN(I),LOCPAR(I),GEO)
         J=LOCDEP(I)
   10 GEO(LOCN,J)=VALUE
      RETURN
      END
      SUBROUTINE THERMO(A,B,C,LINEAR,SYM,WT,VIBS,NVIBS,ESCF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VIBS(*)
      LOGICAL LINEAR
      CHARACTER*80 KEYWRD, KOMENT, TITLE, TMPKEY*80
      COMMON /KEYWRD/ KEYWRD
      COMMON /TITLES/ KOMENT,TITLE
C
C
C   THERMO CALCULATES THE VARIOUS THERMODYNAMIC QUANTITIES FOR A
C   SPECIFIED TEMPERATURE GIVEN THE VIBRATIONAL FREQUENCIES, MOMENTS OF
C   INERTIA, MOLECULAR WEIGHT AND SYMMETRY NUMBER.
C
C   REFERENCE: G.HERZBERG MOLECULAR SPECTRA AND MOLECULAR STRUCTURE
C              VOL 2, CHAP. 5
C
C   ----    TABLE OF SYMMETRY NUMBERS    ----
C
C        C1 CI CS     1      D2 D2D D2H  4       C(INF)V   1
C        C2 C2V C2H   2      D3 D3D D3H  6       D(INF)H   2
C        C3 C3V C3H   3      D4 D4D D4H  8       T TD     12
C        C4 C4V C4H   4      D6 D6D D6H  12      OH       24
C        C6 C6V C6H   6      S6          3
C
C
C   PROGRAM LIMITATIONS:  THE EQUATIONS USED ARE APPROPRIATE TO THE
C   HIGH TEMPERATURE LIMIT AND WILL BEGIN TO BE INADEQUATE AT TEMPERA-
C   TURES BELOW ABOUT 100 K.  SECONDLY THIS PROGRAM IS ONLY APPROPRIATE
C   IN THE CASE OF MOLECULES IN WHICH THERE IS NO FREE ROTATION
C
C
C
C
*******************************************************************
*
*  THE FOLLOWING CONSTANTS ARE NOW DEFINED:
*          PI  = CIRCUMFERENCE TO DIAMETER OF A CIRCLE
*          R   = GAS CONSTANT IN CALORIES/MOLE
*          H   = PLANCK'S CONSTANT IN ERG-SECONDS
*          AK  = BOLTZMANN CONSTANT IN ERG/DEGREE
*          AC  = SPEED OF LIGHT IN CM/SEC
*******************************************************************
      DIMENSION TRANGE(300)
         SAVE                                                           GL0892
      DATA PI /3.14159D0 /
      DATA R/1.98726D0/
      DATA H/6.626D-27/
      DATA AK/1.3807D-16/
      DATA AC/2.99776D+10/
*******************************************************************
      IT1=200
      IT2=400
      ISTEP=10
      TMPKEY=KEYWRD
      I=INDEX(TMPKEY,'THERMO(')
      IF(I.NE.0) THEN
C
C   ERASE ALL TEXT FROM TMPKEY EXCEPT THERMO DATA
C
         TMPKEY(:I)=' '
         TMPKEY(INDEX(TMPKEY,')'):)=' '
         IT1=READA(TMPKEY,I)
         IF(IT1.LT.100) THEN
            WRITE(6,'(//10X,''TEMPERATURE RANGE STARTS TOO LOW,'',
     1'' LOWER BOUND IS RESET TO 30K'')')
            IT1=100
         ENDIF
         I=INDEX(TMPKEY,',')
         IF(I.NE.0) THEN
            TMPKEY(I:I)=' '
            IT2=READA(TMPKEY,I)
            IF(IT2.LT.IT1) THEN
               IT2=IT1+200
               ISTEP=10
               GOTO 10
            ENDIF
            I=INDEX(TMPKEY,',')
            IF(I.NE.0) THEN
               TMPKEY(I:I)=' '
               ISTEP=READA(TMPKEY,I)
               IF(ISTEP.LT.1)ISTEP=1
            ELSE
               ISTEP=(IT2-IT1)/20
               IF(ISTEP.EQ.0)ISTEP=1
               IF(ISTEP.GE.2.AND. ISTEP.LT.5)ISTEP=2
               IF(ISTEP.GE.5.AND. ISTEP.LT.10)ISTEP=5
               IF(ISTEP.GE.10.AND. ISTEP.LT.20)ISTEP=10
               IF(ISTEP.GT.20.AND. ISTEP.LT.50)ISTEP=20
               IF(ISTEP.GT.50.AND. ISTEP.LT.100)ISTEP=50
               IF(ISTEP.GT.100)ISTEP=100
            ENDIF
         ELSE
            IT2=IT1+200
         ENDIF
      ENDIF
   10 CONTINUE
      WRITE(6,'(//,A)')TITLE
      WRITE(6,'(A)')KOMENT
      IF(LINEAR) THEN
         WRITE(6,'(//10X,''MOLECULE IS LINEAR'')')
      ELSE
         WRITE(6,'(//10X,''MOLECULE IS NOT LINEAR'')')
      ENDIF
      WRITE(6,'(/10X,''THERE ARE'',I3,'' GENUINE VIBRATIONS IN THIS '',
     1''SYSTEM'')')NVIBS
      WRITE(6,20)
   20 FORMAT(10X,'THIS THERMODYNAMICS CALCULATION IS LIMITED TO',/
     110X,'MOLECULES WHICH HAVE NO INTERNAL ROTATIONS'//)
      WRITE(6,'(//20X,''CALCULATED THERMODYNAMIC PROPERTIES'')')
      WRITE(6,'(42X,''*'')')
      WRITE(6,'(''   TEMP. (K)   PARTITION FUNCTION   H.O.F.'',
     1''    ENTHALPY   HEAT CAPACITY  ENTROPY'')')
      WRITE(6,'(  ''                                    KCAL/MOL'',
     1''   CAL/MOLE    CAL/K/MOL   CAL/K/MOL'',/)')
      DO 30 I=1,NVIBS
   30 VIBS(I)=ABS(VIBS(I))
      ILIM=1
      DO 40 ITEMP=IT1,IT2,ISTEP
         ILIM=ILIM+1
   40 TRANGE(ILIM)=ITEMP
      TRANGE(1)=298.D0
      DO 80 IR=1,ILIM
         ITEMP=TRANGE(IR)
         T=ITEMP
C   ***   INITIALISE SOME VARIABLES   ***
         C1=H*AC/AK/T
         QV=1.0D0
         HV=0.0D0
         E0=0.0D0
         CPV=0.0D0
         SV1=0.0D0
         SV2=0.0D0
C   ***   CONSTRUCT THE FREQUENCY DEPENDENT PARTS OF PARTITION FUNCTION
         DO 50 I=1,NVIBS
            WI=VIBS(I)
            EWJ=EXP(-WI*C1)
            QV=QV/(1-EWJ)
            HV=HV+WI*EWJ/(1-EWJ)
            E0=E0+WI
            CPV=CPV+WI*WI*EWJ/(1-EWJ)/(1-EWJ)
            SV1=SV1+LOG(1.0D0-EWJ)
   50    SV2=SV2+WI*EWJ/(1-EWJ)
C   ***   FINISH CALCULATION OF VIBRATIONAL PARTS   ***
         HV=HV*R*H*AC/AK
         E0=E0*1.4295D0
         CPV=CPV*R*C1*C1
         SV=SV2*R*C1-R*SV1
C   ***   NOW CALCULATE THE ROTATIONAL PARTS  (FIRST LINEAR MOLECULES
         IF(.NOT.LINEAR) GOTO 60
         QR=1/(C1*A*SYM)
         HR=R*T
         CPR=R
         SR=R*(LOG(T*AK/(H*AC*A*SYM)))+R
         GOTO 70
   60    QR=SQRT(PI/(A*B*C*C1*C1*C1))/SYM
         HR=3.0D0*R*T/2.0D0
         CPR=3.0D0*R/2.0D0
         SR=0.5D0*R*(3.D0*LOG(T*AK/(H*AC))
     1-2.D0*LOG(SYM)+LOG(PI/(A*B*C))+3.D0)
   70    CONTINUE
C   ***   CALCULATE INTERNAL CONTRIBUTIONS   ***
         QINT=QV*QR
         HINT=HV+HR
         CPINT=CPV+CPR
         SINT=SV+SR
C   ***   CONSTRUCT TRANSLATION CONTRIBUTIONS   ***
         QTR=(SQRT(2.D0*PI*WT*T*AK*1.6606D-24)/H)**3
         HTR=5.0D0*R*T/2.0D0
         CPTR=5.0D0*R/2.0D0
         STR=2.2868D0*(5.0D0*LOG10(T)+3.0D0*LOG10(WT))-2.3135D0
C   ***   CONSTRUCT TOTALS   ***
         CPTOT=CPTR+CPINT
         STOT=STR+SINT
         HTOT=HTR+HINT
C   ***   OUTPUT SECTION   ***
         IF(IR.EQ.1)THEN
            H298=HTOT
         ELSE
            WRITE(6,'(/,I7,''  VIB.'',G18.4
     1           ,13X,3F11.5        )')ITEMP,QV,  HV,  CPV,  SV
            WRITE(6,'(7X,''  ROT.'',G13.3
     1           ,16X,3F11.3        )')      QR,  HR,  CPR,  SR
            WRITE(6,'(7X,''  INT.'',G13.3
     1           ,16X,3F11.3        )')      QINT,HINT,CPINT,SINT
            WRITE(6,'(7X,''  TRA.'',G13.3
     1           ,16X,3F11.3)')
     2                                      QTR, HTR, CPTR, STR
            WRITE(6,'(7X,''  TOT.'',13X,F17.3,F11.4,2F11.4)')
     1                     ESCF+(HTOT-H298)/1000.D0,HTOT,CPTOT,STOT
         ENDIF
   80 CONTINUE
      WRITE(6,'(/3X,'' * NOTE: HEATS OF FORMATION ARE RELATIVE TO THE'',
     1/12X,'' ELEMENTS IN THEIR STANDARD STATE AT 298K'')')
      END
C
      SUBROUTINE TIMOUT(NOUT,TIM)
C
C     CONVERT THE TIME FROM SECONDS TO DAYS, HOURS, MINUTES, AND SECONDS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DOUBLE PRECISION MINS, MINPHR
         SAVE                                                           GL0892
C
C
      DATA HRSPD /24.0D0/, MINPHR /60.0D0/
      DATA SECPD /86400.0D0/, SECPHR /3600.0D0/, SECPMI /60.0D0/
C
      DAYS = TIM / SECPD
      IDAYS = INT(DAYS)
      HOURS = (DAYS - FLOAT(IDAYS)) * HRSPD
      IHOURS = INT(HOURS)
      MINS = (HOURS - FLOAT(IHOURS)) * MINPHR
      IMINS = INT(MINS)
      SECS = (MINS - FLOAT(IMINS)) * SECPMI
C
      IF (IDAYS .GT. 1) THEN
         WRITE (NOUT,10) IDAYS,IHOURS,IMINS,SECS
      ELSE IF (IDAYS .EQ. 1) THEN
         WRITE (NOUT,20) IDAYS,IHOURS,IMINS,SECS
      ELSE IF (IHOURS .GT. 0) THEN
         WRITE (NOUT,30) IHOURS,IMINS,SECS
      ELSE IF (IMINS .GT. 0) THEN
         WRITE (NOUT,40) IMINS,SECS
      ELSE
         WRITE (NOUT,50) SECS
      END IF
C
   10 FORMAT (10X,'COMPUTATION TIME = ',I2,1X,'DAYS',2X,I2,1X,'HOURS',
     1        1X,I2,1X,'MINUTES AND',1X,F7.3,1X,'SECONDS')
   20 FORMAT (10X,'COMPUTATION TIME = ',I2,1X,'DAY',2X,I2,1X,'HOURS',
     1        1X,I2,1X,'MINUTES AND',1X,F7.3,1X,'SECONDS')
   30 FORMAT (10X,'COMPUTATION TIME = ',I2,1X,'HOURS',
     1        1X,I2,1X,'MINUTES AND',1X,F7.3,1X,'SECONDS')
   40 FORMAT (10X,'COMPUTATION TIME = ',I2,1X,'MINUTES AND',
     1        1X,F7.3,1X,'SECONDS')
   50 FORMAT (10X,'COMPUTATION TIME = ',F7.3,1X,'SECONDS')
      END
      SUBROUTINE UPDATE(IPARAM, IELMNT, PARAM, MODE,KFN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
************************************************************************
*
*  UPDATE UPDATES THE COMMON BLOCKS WHICH HOLD ALL THE PARAMETERS FOR
*         RUNNING MNDO.
*         IPARAM REFERS TO THE TYPE OF PARAMETER,
*         IELMNT REFERS TO THE ELEMENT,
*         PARAM IS THE VALUE OF THE PARAMETER, AND
*         IF MODE = 1 THEN A COMMON BLOCK IS UPDATED,
*         IF MODE = 2 THEN A DATUM IS EXTRACTED FROM THE COMMON BLOCK.
*
************************************************************************
      COMMON /MNDO/  USSM(107), UPPM(107), UDDM(107), ZSM(107),ZPM(107),
     1ZDM(107), BETASM(107), BETAPM(107), BETADM(107), ALPM(107),
     2EISOLM(107), DDM(107), QQM(107), AMM(107), ADM(107), AQM(107)
     3,GSSM(107),GSPM(107),GPPM(107),GP2M(107),HSPM(107), POLVOM(107)
      COMMON /EXPONT/ ZS(107),ZP(107),ZD(107)
     1       /NATORB/ NATORB(107)
     2       /BETAS / BETAS(107),BETAP(107),BETAD(107)
     3       /VSIPS / VS(107),VP(107),VD(107)
     4       /ONELEC/ USS(107),UPP(107),UDD(107)
     5       /MULTIP/ DD(107),QQ(107),AM(107),AD(107),AQ(107)
     6       /TWOELE/ GSS(107),GSP(107),GPP(107),GP2(107),HSP(107)
     7                ,GSD(107),GPD(107),GDD(107)
     8       /ALPHA / ALP(107)
     9       /IDEAS / GUESS1(107,10), GUESS2(107,10), GUESS3(107,10)
      COMMON /GAUSS / FN1(107),FN2(107)
         SAVE                                                           GL0892
      GOTO
     1(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,
     2190,200,210,220,230,240,250),IPARAM
   10 USS (IELMNT)=PARAM
      USSM(IELMNT)=PARAM
      RETURN
   20 UPP (IELMNT)=PARAM
      UPPM(IELMNT)=PARAM
      RETURN
   30 UDD (IELMNT)=PARAM
      UDDM(IELMNT)=PARAM
      RETURN
   40 ZS (IELMNT)=PARAM
      ZSM(IELMNT)=PARAM
      RETURN
   50 ZP (IELMNT)=PARAM
      ZPM(IELMNT)=PARAM
      RETURN
   60 ZD (IELMNT)=PARAM
      ZDM(IELMNT)=PARAM
      RETURN
   70 BETAS (IELMNT)=PARAM
      BETASM(IELMNT)=PARAM
      RETURN
   80 BETAP (IELMNT)=PARAM
      BETAPM(IELMNT)=PARAM
      RETURN
   90 BETAD (IELMNT)=PARAM
      BETADM(IELMNT)=PARAM
      RETURN
  100 GSS (IELMNT)=PARAM
      GSSM(IELMNT)=PARAM
      RETURN
  110 GSP (IELMNT)=PARAM
      GSPM(IELMNT)=PARAM
      RETURN
  120 GPP (IELMNT)=PARAM
      GPPM(IELMNT)=PARAM
      RETURN
  130 GP2 (IELMNT)=PARAM
      GP2M(IELMNT)=PARAM
      RETURN
  140 HSP (IELMNT)=PARAM
      HSPM(IELMNT)=PARAM
      RETURN
  150 RETURN
  160 RETURN
  170 RETURN
  180 ALP (IELMNT)=PARAM
      ALPM(IELMNT)=PARAM
      RETURN
  190 RETURN
  200 RETURN
  210 RETURN
  220 GUESS1(IELMNT,KFN)=PARAM
      RETURN
  230 GUESS2(IELMNT,KFN)=PARAM
      RETURN
  240 GUESS3(IELMNT,KFN)=PARAM
      RETURN
  250 NATORB(IELMNT)=PARAM
      I=INT(PARAM+0.5)
      IF(I.NE.9.AND.I.NE.4.AND.I.NE.1)THEN
         WRITE(6,'(///10X,'' UNACCEPTABLE VALUE FOR NO. OF ORBITALS'',
     1'' ON ATOM'')')
         STOP
      ENDIF
      END
      SUBROUTINE VECPRT (A,NUMM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"
C
C
      DIMENSION  A(*)
c Common MOLKST splitted in MOLKSI and MOLKSR    Ivan Rossi 0394   &8)
      COMMON /MOLKSI/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),
     1                NMIDLE(NUMATM),NLAST(NUMATM), NORBS,
     2                NELECS,NALPHA,NBETA,NCLOSE,NOPEN
      COMMON /MOLKSR/ FRACT
      COMMON /ELEMTS/ ELEMNT(107)
C**********************************************************************
C
C  VECPRT PRINTS A LOWER-HALF TRIANGLE OF A SQUARE MATRIX, THE
C         LOWER-HALF TRIANGLE BEING STORED IN PACKED FORM IN THE
C         ARRAY "A"
C
C ON INPUT:
C      A      = ARRAY TO BE PRINTED
C      NUMM   = SIZE OF ARRAY TO BE PRINTED
C(REF) NUMAT  = NUMBER OF ATOMS IN THE MOLECULE (THIS IS NEEDED TO
C               DECIDE IF AN ATOMIC ARRAY OR ATOMIC ORBITAL ARRAY IS
C               TO BE PRINTED
C(REF) NAT    = LIST OF ATOMIC NUMBERS
C(REF) NFIRST = LIST OF ORBITAL COUNTERS
C(REF) NLAST  = LIST OF ORBITAL COUNTERS
C
C  NONE OF THE ARGUMENTS ARE ALTERED BY THE CALL OF VECPRT
C
C*********************************************************************
      DIMENSION NATOM(MAXORB)
      CHARACTER * 6 LINE(21)
      CHARACTER*2 ELEMNT,ATORBS(9), ITEXT(MAXORB), JTEXT(MAXORB)
         SAVE                                                           GL0892
      DATA ATORBS/' S','PX','PY','PZ','X2','XZ','Z2','YZ','XY'/
      IF(NUMAT.NE.0.AND.NUMAT.EQ.NUMM) THEN
C
C    PRINT OVER ATOM COUNT
C
         DO 10 I=1,NUMAT
            ITEXT(I)='  '
            JTEXT(I)=ELEMNT(NAT(I))
            NATOM(I)=I
   10    CONTINUE
      ELSE
         IF (NUMAT.NE.0.AND.NLAST(NUMAT) .EQ. NUMM) THEN
            DO 30 I=1,NUMAT
               JLO=NFIRST(I)
               JHI=NLAST(I)
               L=NAT(I)
               K=0
               DO 20 J=JLO,JHI
                  K=K+1
                  ITEXT(J)=ATORBS(K)
                  JTEXT(J)=ELEMNT(L)
                  NATOM(J)=I
   20          CONTINUE
   30       CONTINUE
         ELSE
            NUMB=ABS(NUMM)
            DO 40 I=1,NUMB
               ITEXT(I) = '  '
               JTEXT(I) = '  '
   40       NATOM(I)=I
         ENDIF
      ENDIF
      NUMB=ABS(NUMM)
      DO 50 I=1,21
   50 LINE(I)='------'
      LIMIT=(NUMB*(NUMB+1))/2
      KK=8
      NA=1
   60 LL=0
      M=MIN0((NUMB+1-NA),6)
      MA=2*M+1
      M=NA+M-1
      WRITE(6,100)(ITEXT(I),JTEXT(I),NATOM(I),I=NA,M)
      WRITE (6,110) (LINE(K),K=1,MA)
      DO 80 I=NA,NUMB
         LL=LL+1
         K=(I*(I-1))/2
         L=MIN0((K+M),(K+I))
         K=K+NA
         IF ((KK+LL).LE.50) GO TO 70
         WRITE (6,120)
         WRITE (6,100) (ITEXT(N),JTEXT(N),NATOM(N),N=NA,M)
         WRITE (6,110) (LINE(N),N=1,MA)
         KK=4
         LL=0
   70    WRITE (6,130) ITEXT(I),JTEXT(I),NATOM(I),(A(N),N=K,L)
   80 CONTINUE
      IF (L.GE.LIMIT) GO TO 90
      KK=KK+LL+4
      NA=M+1
      IF ((KK+NUMB+1-NA).LE.50) GO TO 60
      KK=4
      WRITE (6,120)
      GO TO 60
   90 RETURN
C
  100 FORMAT (1H0/13X,10(1X,A2,1X,A2,I3,2X))
  110 FORMAT (1H ,21A6)
  120 FORMAT (1H1)
  130 FORMAT (1H ,A2,1X,A2,I5,10F11.6)
C
      END
      LOGICAL FUNCTION WORD(KEYWRD,TESTWD)
      CHARACTER KEYWRD*80, TESTWD*(*)
         SAVE                                                           GL0892
      WORD=.FALSE.
   10 J=INDEX(KEYWRD,TESTWD)
      IF(J.NE.0)THEN
   20    IF(KEYWRD(J:J).NE.' ')GOTO 30
         J=J+1
         GOTO 20
   30    WORD=.TRUE.
         DO 60 K=J,80
            IF(KEYWRD(K:K).EQ.'='.OR.KEYWRD(K:K).EQ.' ') THEN
C
C     CHECK FOR ATTACHED '=' SIGN
C
               J=K
               IF(KEYWRD(J:J).EQ.'=')GOTO 50
C
C     CHECK FOR SEPARATED '=' SIGN
C
               DO 40 J=K+1,80
                  IF(KEYWRD(J:J).EQ.'=') GOTO 50
   40          IF(KEYWRD(J:J).NE.' ')GOTO 10
C
C    THERE IS NO '=' SIGN ASSOCIATED WITH THIS KEYWORD
C
               GOTO 10
   50          KEYWRD(J:J)=' '
C
C   THERE MUST BE A NUMBER AFTER THE '=' SIGN, SOMEWHERE
C
               GOTO 20
            ENDIF
   60    KEYWRD(K:K)=' '
      ENDIF
      RETURN
      END
      SUBROUTINE XYZGEO(XYZ,NUMAT,NA,NB,NC,DEGREE,GEO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XYZ(3,*), NA(*), NB(*), NC(*), GEO(3,*)
         SAVE                                                           GL0892
***********************************************************************
*
*   XYZGEO CONVERTS COORDINATES FROM CARTESIAN TO INTERNAL.
*
*     ON INPUT XYZ  = ARRAY OF CARTESIAN COORDINATES
*              NUMAT= NUMBER OF ATOMS
*              NA   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
*                     BY DISTANCE
*              NB   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
*                     BY ANGLE
*              NC   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
*                     BY DIHEDRAL
*
*    ON OUTPUT GEO  = INTERNAL COORDINATES IN ANGSTROMS, RADIANS,
*                     AND RADIANS
*
***********************************************************************
      DO 30 I=2,NUMAT
         J=NA(I)
         K=NB(I)
         L=NC(I)
         IF(I.LT.3) GOTO 30
         II=I
         CALL BANGLE(XYZ,II,J,K,GEO(2,I))
         GEO(2,I)=GEO(2,I)*DEGREE
         IF(I.LT.4) GOTO 30
C
C   MAKE SURE DIHEDRAL IS MEANINGLFUL
C
         CALL BANGLE(XYZ,J,K,L,ANGL)
         TOL=0.2617994D0
         IF(ANGL.GT.3.1415926D0-TOL.OR.ANGL.LT.TOL)THEN
C
C  ANGLE IS UNSATISFACTORY, LET'S SEARCH FOR ANOTHER ATOM FOR
C  DEFINING THE DIHEDRAL.
   10       SUM=100.D0
            DO 20 I1=1,II-1
               R=(XYZ(1,I1)-XYZ(1,K))**2+
     1          (XYZ(2,I1)-XYZ(2,K))**2+
     2          (XYZ(3,I1)-XYZ(3,K))**2
               IF(R.LT.SUM.AND.I1.NE.J.AND.I1.NE.K) THEN
                  CALL BANGLE(XYZ,J,K,I1,ANGL)
                  IF(ANGL.LT.3.1415926D0-TOL.AND.ANGL.GT.TOL)THEN
                     SUM=R
                     L=I1
                     NC(II)=L
                  ENDIF
               ENDIF
   20       CONTINUE
            IF(SUM.GT.99.D0.AND.TOL.GT.0.1D0)THEN
C
C ANYTHING WITHIN 5 DEGREES?
C
               TOL=0.087266D0
               GOTO 10
            ENDIF
         ENDIF
         CALL DIHED(XYZ,II,J,K,L,GEO(3,I))
         GEO(3,I)=GEO(3,I)*DEGREE
   30 GEO(1,I)= SQRT((XYZ(1,I)-XYZ(1,J))**2+
     1                   (XYZ(2,I)-XYZ(2,J))**2+
     2                   (XYZ(3,I)-XYZ(3,J))**2)
      GEO(1,1)=0.D0
      GEO(2,1)=0.D0
      GEO(3,1)=0.D0
      GEO(2,2)=0.D0
      GEO(3,2)=0.D0
      GEO(3,3)=0.D0
      RETURN
      END
      SUBROUTINE XYZINT(XYZ,NUMAT,NA,NB,NC,DEGREE,GEO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XYZ(3,*), NA(*), NB(*), NC(*), GEO(3,*)
***********************************************************************
*
* XYZINT WORKS OUT THE INTERNAL COORDINATES OF A MOLECULE.
*        THE "RULES" FOR THE CONNECTIVITY ARE AS FOLLOWS:
*        ATOM I IS DEFINED AS BEING AT A DISTANCE FROM THE NEAREST
*        ATOM J, ATOM J ALREADY HAVING BEEN DEFINED.
*        ATOM I MAKES AN ANGLE WITH ATOM J AND THE ATOM K, WHICH HAS
*        ALREADY BEEN DEFINED, AND IS THE NEAREST ATOM TO J
*        ATOM I MAKES A DIHEDRAL ANGLE WITH ATOMS J, K, AND L. L HAVING
*        BEEN DEFINED AND IS THE NEAREST ATOM TO K, AND J, K AND L
*        HAVE A CONTAINED ANGLE IN THE RANGE 15 TO 165 DEGREES,
*        IF POSSIBLE.
*
*        IF(NA(2).EQ.1 THEN THE ORIGINAL CONNECTIVITY IS USED.
*
*        NOTE THAT GEO AND XYZ MUST NOT BE THE SAME IN THE CALL.
*
*   ON INPUT XYZ    = CARTESIAN ARRAY OF NUMAT ATOMS
*            DEGREE = 1 IF ANGLES ARE TO BE IN RADIANS
*            DEGREE = 57.29578 IF ANGLES ARE TO BE IN DEGREES
*
***********************************************************************
      COMMON /GEOOK/ IGEOOK
      COMMON /NUMCAL/ NUMCAL                                            IR0295
         SAVE                                                           GL0892
      DATA ICALCN/0/                                                    IR0295
      IGEOOK=99
      NAI1=0
      NAI2=0
      IF((ICALCN.EQ.NUMCAL).AND.NA(2).EQ.-1)THEN                        IR0295
         NA(2)=1
         DO 10 I=2,NUMAT
            J=NA(I)
            IF(I.GT.3)CALL DIHED(XYZ,I,J,NB(I),NC(I),GEO(3,I))
            IF(I.GT.2)CALL BANGLE(XYZ,I,J,NB(I),GEO(2,I))
            GEO(1,I)=SQRT((XYZ(1,I)-XYZ(1,J))**2+
     1          (XYZ(2,I)-XYZ(2,J))**2+
     2          (XYZ(3,I)-XYZ(3,J))**2)
   10    CONTINUE
      ELSE
         IF(NA(2).EQ.-1) ICALCN=NUMCAL                                   IR0295
         DO 30 I=1,NUMAT
            NA(I)=2
            NB(I)=3
            NC(I)=4
            IM1=I-1
            IF(IM1.EQ.0)GOTO 30
            SUM=1.D30
            DO 20 J=1,IM1
               R=(XYZ(1,I)-XYZ(1,J))**2+
     1          (XYZ(2,I)-XYZ(2,J))**2+
     2          (XYZ(3,I)-XYZ(3,J))**2
               IF(R.LT.SUM.AND.NA(J).NE.J.AND.NB(J).NE.J) THEN
                  SUM=R
                  K=J
               ENDIF
   20       CONTINUE
C
C   ATOM I IS NEAREST TO ATOM K
C
            NA(I)=K
            IF(I.GT.2)NB(I)=NA(K)
            IF(I.GT.3)NC(I)=NB(K)
C
C   FIND ANY ATOM TO RELATE TO NA(I)
C
   30    CONTINUE
      ENDIF
      NA(1)=0
      NB(1)=0
      NC(1)=0
      NB(2)=0
      NC(2)=0
      NC(3)=0
      CALL XYZGEO(XYZ,NUMAT,NA,NB,NC,DEGREE,GEO)
      RETURN
      END
