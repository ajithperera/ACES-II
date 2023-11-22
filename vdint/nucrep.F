      SUBROUTINE NUCREP(ENERNN,GRADNN,HESSNN,NCOORD)
C
C     This subroutine calculates the nuclear repulsion contributions
C     to the molecular gradient and Hessian.
C
CEND
C
C     ADAPTED TO THE ACES II ENVIRONMENT OCT/90/JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "mxcent.par"
#include "baslims.par"
      PARAMETER (AZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0,
     *           THREE = 3.0D0)
      LOGICAL SECOND
      LOGICAL SKIP
      LOGICAL DCORD, DCORGD, NOORBT, DOPERT
      CHARACTER NAMEX*6
C
      DIMENSION GRADNN(NCOORD),HESSNN(NCOORD,NCOORD)
C
      COMMON /CBINUC/ JPRINT, SKIP, MAXDIF
      COMMON /NUCLEIi/ NOORBT(MXCENT),
     &                NUCIND, NUCDEP, NUCPRE(MXCENT), NUCNUM(MXCENT,8),
     &                NUCDEG(MXCENT), ISTBNU(MXCENT), NDCORD,
     &                NDCOOR(MXCOOR), NTRACO, NROTCO, ITRACO(3),
     &                IROTCO(3),
     &                NATOMS, NFLOAT,
     &                IPTGDV(3*MXCENT),
     &                NGDVEC(8), IGDVEC(8)
      COMMON /NUCLEI/ CHARGE(MXCENT), CORD(MXCENT,3),
     &                DCORD(MXCENT,3),DCORGD(MXCENT,3),
     &                DOPERT(0:3*MXCENT)
      COMMON /NUCLEC/ NAMEX(MXCOOR)
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
      COMMON/IPRI/IPRINT
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
      IF (SKIP) RETURN
C
C  INITIALIZE 
C
      SECOND = MAXDIF .EQ. 2
      CALL ZERO(GRADNN(1),NCOORD)
      IF (SECOND) CALL ZERO(HESSNN,NCOORD*NCOORD)
C
C  CALCULATE GRADIENT AND HESSIAN 
C
C     Run over symmetry-independent nuclei A
C
      DO 200 NCENTA = 1, NUCIND
         CHARGA = CHARGE(NCENTA)
         IF (ABS(CHARGA) .GT. AZERO) THEN
            NAX    = 3*NCENTA - 2
            NAY    = 3*NCENTA - 1
            NAZ    = 3*NCENTA
            CORDXA = CORD(NCENTA,1)
            CORDYA = CORD(NCENTA,2)
            CORDZA = CORD(NCENTA,3)
            MULA   = ISTBNU(NCENTA)
C
C  Run over symmetry-independent nuclei B
C
            DO 300 NCENTB =  1, NCENTA
               CHARGB = CHARGE(NCENTB)
               IF (ABS(CHARGB) .GT. AZERO) THEN
                  NBX    = 3*NCENTB - 2
                  NBY    = 3*NCENTB - 1
                  NBZ    = 3*NCENTB
                  CORBX0 = CORD(NCENTB,1)
                  CORBY0 = CORD(NCENTB,2)
                  CORBZ0 = CORD(NCENTB,3)
                  MULB   = ISTBNU(NCENTB)
C
                  MAB    = IBTOR (MULA,MULB)
                  KAB    = IBTAND(MULA,MULB)
                  HKAB   = FMULT(KAB)
                  CROSS = ONE
                  IF (NCENTA .EQ. NCENTB) THEN
                     HKAB = HALF*HKAB
                     CROSS = TWO
                  END IF
C
C  Run over symmetry independent charge distributions
C
                  DO 400 ISYMOP = 0, MAXLOT
                  IF (IBTAND(ISYMOP,MAB) .EQ. 0) THEN
                     ICENTA = NUCNUM(NCENTA,1)
                     ICENTB = NUCNUM(NCENTB,ISYMOP+1)
                  IF (ICENTA .EQ. ICENTB) GO TO 400
C
                     SBX = PT(IBTAND(ISYTYP(1),ISYMOP))
                     SBY = PT(IBTAND(ISYTYP(2),ISYMOP))
                     SBZ = PT(IBTAND(ISYTYP(3),ISYMOP))
                     XAB = CORDXA - SBX*CORBX0
                     YAB = CORDYA - SBY*CORBY0
                     ZAB = CORDZA - SBZ*CORBZ0
C
                     XAB2   = XAB*XAB
                     YAB2   = YAB*YAB
                     ZAB2   = ZAB*ZAB
                     RAB2   = XAB2 + YAB2 + ZAB2
                     RAB1   = SQRT(RAB2)
                     ZZR3IN = - HKAB*CHARGA*CHARGB/(RAB1*RAB2)
                     VNUCX  = XAB*ZZR3IN
                     VNUCY  = YAB*ZZR3IN
                     VNUCZ  = ZAB*ZZR3IN
C
C   Gradient 
C
C   Totally symmetric contributions only
C
                     IAX  = IPTCNT(NAX,0)
                     IAY  = IPTCNT(NAY,0)
                     IAZ  = IPTCNT(NAZ,0)
                     IBX  = IPTCNT(NBX,0)
                     IBY  = IPTCNT(NBY,0)
                     IBZ  = IPTCNT(NBZ,0)
                     IF (IAX.NE.0) GRADNN(IAX) = GRADNN(IAX) + VNUCX
                     IF (IAY.NE.0) GRADNN(IAY) = GRADNN(IAY) + VNUCY
                     IF (IAZ.NE.0) GRADNN(IAZ) = GRADNN(IAZ) + VNUCZ
                     IF (IBX.NE.0) GRADNN(IBX) = GRADNN(IBX) - SBX*VNUCX
                     IF (IBY.NE.0) GRADNN(IBY) = GRADNN(IBY) - SBY*VNUCY
                     IF (IBZ.NE.0) GRADNN(IBZ) = GRADNN(IBZ) - SBZ*VNUCZ
C
C  Hessian 
C
                     IF (SECOND) THEN
                        ZZR5IN = ZZR3IN/RAB2
                        VNUCXX = (XAB2 + XAB2 - YAB2 - ZAB2)*ZZR5IN
                        VNUCXY = THREE*XAB*YAB*ZZR5IN
                        VNUCXZ = THREE*XAB*ZAB*ZZR5IN
                        VNUCYY = (YAB2 + YAB2 - ZAB2 - XAB2)*ZZR5IN
                        VNUCYZ = THREE*YAB*ZAB*ZZR5IN
                        VNUCZZ = (ZAB2 + ZAB2 - XAB2 - YAB2)*ZZR5IN
C
C  Loop over irreps for differentiation operator
C
                        DO 500 IREPD = 0, MAXLOP
                           CHI = PT(IBTAND(ISYMOP,IREPD))
                           IAX  = IPTCNT(NAX,IREPD)
                           IAY  = IPTCNT(NAY,IREPD)
                           IAZ  = IPTCNT(NAZ,IREPD)
                           IBX  = IPTCNT(NBX,IREPD)
                           IBY  = IPTCNT(NBY,IREPD)
                           IBZ  = IPTCNT(NBZ,IREPD)
                           IF (IAX*IBX.NE.0) THEN
                              HESSNN(IAX,IBX) = HESSNN(IAX,IBX)
     *                                        + CROSS*CHI*SBX*VNUCXX
                           END IF
                           IF (IAX*IBY.NE.0) THEN
                              HESSNN(IAX,IBY) = HESSNN(IAX,IBY)
     *                                        + CHI*SBY*VNUCXY
                           END IF
                           IF (IAX*IBZ.NE.0) THEN
                              HESSNN(IAX,IBZ) = HESSNN(IAX,IBZ)
     *                                        + CHI*SBZ*VNUCXZ
                           END IF
                           IF (IAY*IBX.NE.0) THEN
                              HESSNN(IAY,IBX) = HESSNN(IAY,IBX)
     *                                        + CHI*SBX*VNUCXY
                           END IF
                           IF (IAY*IBY.NE.0) THEN
                              HESSNN(IAY,IBY) = HESSNN(IAY,IBY)
     *                                        + CROSS*CHI*SBY*VNUCYY
                           END IF
                           IF (IAY*IBZ.NE.0) THEN
                              HESSNN(IAY,IBZ) = HESSNN(IAY,IBZ)
     *                                        + CHI*SBZ*VNUCYZ
                           END IF
                           IF (IAZ*IBX.NE.0) THEN
                              HESSNN(IAZ,IBX) = HESSNN(IAZ,IBX)
     *                                        + CHI*SBX*VNUCXZ
                           END IF
                           IF (IAZ*IBY.NE.0) THEN
                              HESSNN(IAZ,IBY) = HESSNN(IAZ,IBY)
     *                                        + CHI*SBY*VNUCYZ
                           END IF
                           IF (IAZ*IBZ.NE.0) THEN
                              HESSNN(IAZ,IBZ) = HESSNN(IAZ,IBZ)
     *                                        + CROSS*CHI*SBZ*VNUCZZ
                           END IF
                           IF (IAX*IAX.NE.0) THEN
                              HESSNN(IAX,IAX) = HESSNN(IAX,IAX) - VNUCXX
                           END IF
                           IF (IAY*IAX.NE.0) THEN
                              HESSNN(IAY,IAX) = HESSNN(IAY,IAX) - VNUCXY
                           END IF
                           IF (IAY*IAY.NE.0) THEN
                              HESSNN(IAY,IAY) = HESSNN(IAY,IAY) - VNUCYY
                           END IF
                           IF (IAZ*IAX.NE.0) THEN
                              HESSNN(IAZ,IAX) = HESSNN(IAZ,IAX) - VNUCXZ
                           END IF
                           IF (IAZ*IAY.NE.0) THEN
                              HESSNN(IAZ,IAY) = HESSNN(IAZ,IAY) - VNUCYZ
                           END IF
                           IF (IAZ*IAZ.NE.0) THEN
                              HESSNN(IAZ,IAZ) = HESSNN(IAZ,IAZ) - VNUCZZ
                           END IF
                           IF (IBX*IBX.NE.0) THEN
                              HESSNN(IBX,IBX) = HESSNN(IBX,IBX) - VNUCXX
                           END IF
                           IF (IBY*IBX.NE.0) THEN
                              HESSNN(IBY,IBX) = HESSNN(IBY,IBX)
     *                                        - SBX*SBY*VNUCXY
                           END IF
                           IF (IBY*IBY.NE.0) THEN
                              HESSNN(IBY,IBY) = HESSNN(IBY,IBY) - VNUCYY
                           END IF
                           IF (IBZ*IBX.NE.0) THEN
                              HESSNN(IBZ,IBX) = HESSNN(IBZ,IBX)
     *                                        - SBX*SBZ*VNUCXZ
                           END IF
                           IF (IBZ*IBY.NE.0) THEN
                              HESSNN(IBZ,IBY) = HESSNN(IBZ,IBY)
     *                                        - SBY*SBZ*VNUCYZ
                           END IF
                           IF (IBZ*IBZ.NE.0) THEN
                              HESSNN(IBZ,IBZ) = HESSNN(IBZ,IBZ) - VNUCZZ
                           END IF
  500                   CONTINUE
                     END IF
                  END IF
  400             CONTINUE
               END IF
  300       CONTINUE
         END IF
  200 CONTINUE
C
C  PRINT SECTION 
C
      IF (IPRINT .GT. 1) THEN
       CALL HEADER('Nuclear repulsion gradient',-1)
       CALL PRIGRD(GRADNN,0)
      END IF
      IF (SECOND) THEN
C
C SYMMETRIZE HESSIAN
C
       DO 800 I=1,NCOORD
        DO 800 J=1,I-1
         Z=HESSNN(J,I)+HESSNN(I,J)
         HESSNN(J,I)=Z
         HESSNN(I,J)=Z
800    CONTINUE
       IF (IPRINT .GT. 1) THEN
        CALL HEADER('Nuclear repulsion Hessian',-1)
        CALL PRIHES(HESSNN,NCOORD,'CENTERS',0)
       END IF
      END IF
      RETURN
      END
