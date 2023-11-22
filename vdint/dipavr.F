      SUBROUTINE DIPAVR(KHKTA,KHKTB,IDENA,IDENB,NCENTA,NCENTB,
     *                  ISYMOP,SBX,SBY,SBZ,LDIAG,ONECEN,DIFDIP,
     *                  WORK1,DIPME,DDIPE,NCOORD)
C
C     CALCULATES THE VARIOUS EXPECTATION VALUES FOR DIPOLE MOMENT 
C     AND DIPOLE MOMENT DERIVATIVES
C
CEND
C     tuh 1985
C     symmetry 081288 tuh
C     error fixed 28/10/89 jg
C
C  ADAPTED TO THE ACES II ENVIRONMENT OCT/90/JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "mxcent.par"
#include "baslims.par"
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN*
     *           MXCONT*MXCONT)
      PARAMETER (D0 = 0.00 D00)
      LOGICAL DIFDIP, ONECEN, LDIAG
      DIMENSION WORK1(1),DIPME(3),DDIPE(3,NCOORD)
      COMMON/DF/ DSHELL(MXAQNS),FSHELL(MXAQNS)
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
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
      COMMON/DENPOIN/KDEN,KFOC
      COMMON /GENCON/ NRCA,NRCB,CONTA(MXCONT),CONTB(MXCONT)
C
      DATA AZERO /0.0D0/
C
      IBTAND(I,J) = AND(I,J)
      ITRI(I,J) = MAX(I,J)*(MAX(I,J) - 1)/2 + MIN(I,J)
C
      NAX = 3*NCENTA - 2
      NAY = 3*NCENTA - 1
      NAZ = 3*NCENTA
      NBX = 3*NCENTB - 2
      NBY = 3*NCENTB - 1
      NBZ = 3*NCENTB
      IX  = IPTAX(1)
      IY  = IPTAX(2)
      IZ  = IPTAX(3)
C
C     Density matrix
C
      IOFF=0
      ISKIP=NRCA*NRCB
      IF(LDIAG) ISKIP=NRCA*(NRCA+1)/2
      MAXCMP = ISKIP*KHKTA*KHKTB
C
C     Expectation value of undifferentiated elements
C
      XAVR0=AZERO
      YAVR0=AZERO
      ZAVR0=AZERO
#ifdef _DEBUG_LVLM
      Write(6,*) ISYTYP(1), ISYTYP(2), ISYTYP(3)
      Write(6,*) "X,Y,Z Dipole elements"
      Write(6, "(3(1x,F20.10))") (Work1(ID000x+i),i=1,MAXCMP)
      Write(6,*)
      Write(6, "(3(1x,F20.10))") (Work1(ID000y+i),i=1,MAXCMP)
      Write(6,*)
      Write(6, "(3(1x,F20.10))") (Work1(ID000z+i),i=1,MAXCMP)
      Write(6,*)
      Write(6,*) "density"
      Write(6, "(3(1x,F20.10))") (DSHELL(i),i=1,MAXCMP)
#endif 

      SAVR0=SDOT(MAXCMP,DSHELL,1,WORK1(ID0000+1),1)
      IF(ISYTYP(1).EQ.0)
     *    XAVR0=SDOT(MAXCMP,DSHELL,1,WORK1(ID000X+1),1)
      IF(ISYTYP(2).EQ.0)
     *     YAVR0=SDOT(MAXCMP,DSHELL,1,WORK1(ID000Y+1),1)
      IF(ISYTYP(3).EQ.0)
     *     ZAVR0=SDOT(MAXCMP,DSHELL,1,WORK1(ID000Z+1),1) 
      DIPME(1) = DIPME(1) - XAVR0
      DIPME(2) = DIPME(2) - YAVR0
      DIPME(3) = DIPME(3) - ZAVR0
      IF (DIFDIP) THEN
       IF (ONECEN) THEN
        DO 300 IREP = 0, MAXLOP
         IAX = IPTCNT(NAX,IREP)
         IAY = IPTCNT(NAY,IREP)
         IAZ = IPTCNT(NAZ,IREP)
         IF (ISYTYP(1) .EQ. IREP .AND. IAX .GT. 0) THEN
          DDIPE(IX,IAX) = DDIPE(IX,IAX) - SAVR0
         END IF
         IF (ISYTYP(2) .EQ. IREP .AND. IAY .GT. 0) THEN
          DDIPE(IY,IAY) = DDIPE(IY,IAY) - SAVR0
         END IF
         IF (ISYTYP(3) .EQ. IREP .AND. IAZ .GT. 0) THEN
          DDIPE(IZ,IAZ) = DDIPE(IZ,IAZ) - SAVR0
         END IF
  300   CONTINUE
       ELSE
        XAVRX=-SDOT(MAXCMP,DSHELL,1,WORK1(ID00XX+1),1)
        XAVRY=-SDOT(MAXCMP,DSHELL,1,WORK1(ID00YX+1),1)
        XAVRZ=-SDOT(MAXCMP,DSHELL,1,WORK1(ID00ZX+1),1)
        YAVRX=-SDOT(MAXCMP,DSHELL,1,WORK1(ID00XY+1),1)
        YAVRY=-SDOT(MAXCMP,DSHELL,1,WORK1(ID00YY+1),1)
        YAVRZ=-SDOT(MAXCMP,DSHELL,1,WORK1(ID00ZY+1),1)
        ZAVRX=-SDOT(MAXCMP,DSHELL,1,WORK1(ID00XZ+1),1)
        ZAVRY=-SDOT(MAXCMP,DSHELL,1,WORK1(ID00YZ+1),1)
        ZAVRZ=-SDOT(MAXCMP,DSHELL,1,WORK1(ID00ZZ+1),1)
        DO 500 IREP = 0, MAXLOP
         CHI  = PT(IBTAND(ISYMOP,IREP))
         CSBX = CHI*SBX
         CSBY = CHI*SBY
         CSBZ = CHI*SBZ
         IAX = IPTCNT(NAX,IREP)
         IAY = IPTCNT(NAY,IREP)
         IAZ = IPTCNT(NAZ,IREP)
         IBX = IPTCNT(NBX,IREP)
         IBY = IPTCNT(NBY,IREP)
         IBZ = IPTCNT(NBZ,IREP)
         IF (ISYTYP(1) .EQ. IREP) THEN
          IF (IAX.GT.0)
     *    DDIPE(IX,IAX) = DDIPE(IX,IAX) + XAVRX
          IF (IAY.GT.0)
     *    DDIPE(IX,IAY) = DDIPE(IX,IAY) + XAVRY
          IF (IAZ.GT.0)
     *    DDIPE(IX,IAZ) = DDIPE(IX,IAZ) + XAVRZ
          IF (IBX.GT.0)
     *    DDIPE(IX,IBX) = DDIPE(IX,IBX) - CSBX*XAVRX - SAVR0
          IF (IBY.GT.0)
     *    DDIPE(IX,IBY) = DDIPE(IX,IBY) - CSBY*XAVRY
          IF (IBZ.GT.0)
     *    DDIPE(IX,IBZ) = DDIPE(IX,IBZ) - CSBZ*XAVRZ
         END IF
         IF (ISYTYP(2) .EQ. IREP) THEN
          IF (IAX.GT.0)
     *    DDIPE(IY,IAX) = DDIPE(IY,IAX) + YAVRX
          IF (IAY.GT.0)
     *    DDIPE(IY,IAY) = DDIPE(IY,IAY) + YAVRY
          IF (IAZ.GT.0)
     *    DDIPE(IY,IAZ) = DDIPE(IY,IAZ) + YAVRZ
          IF (IBX.GT.0)
     *    DDIPE(IY,IBX) = DDIPE(IY,IBX) - CSBX*YAVRX
          IF (IBY.GT.0)
     *    DDIPE(IY,IBY) = DDIPE(IY,IBY) - CSBY*YAVRY - SAVR0
          IF (IBZ.GT.0)
     *    DDIPE(IY,IBZ) = DDIPE(IY,IBZ) - CSBZ*YAVRZ
         END IF
         IF (ISYTYP(3) .EQ. IREP) THEN
          IF (IAX.GT.0)
     *    DDIPE(IZ,IAX) = DDIPE(IZ,IAX) + ZAVRX
          IF (IAY.GT.0)
     *    DDIPE(IZ,IAY) = DDIPE(IZ,IAY) + ZAVRY
          IF (IAZ.GT.0)
     *    DDIPE(IZ,IAZ) = DDIPE(IZ,IAZ) + ZAVRZ
          IF (IBX.GT.0)
     *    DDIPE(IZ,IBX) = DDIPE(IZ,IBX) - CSBX*ZAVRX
          IF (IBY.GT.0)
     *    DDIPE(IZ,IBY) = DDIPE(IZ,IBY) - CSBY*ZAVRY
          IF (IBZ.GT.0)
     *    DDIPE(IZ,IBZ) = DDIPE(IZ,IBZ) - CSBZ*ZAVRZ - SAVR0
         END IF
  500   CONTINUE
       END IF
      END IF
      RETURN
      END
