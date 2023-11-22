      SUBROUTINE AVEKFS(KHKTA,KHKTB,IDENA,IDENB,ICENTA,ICENTB,
     &                  NCENTA,NCENTB,SBX,SBY,SBZ,
     &                  ISYMOP,SECOND,LDIAG,WORK1,ENERKE,
     &                  GRADKE,GRADFS,HESSKE,HESFS2,NCOORD)
C
CEND
C
C     ADAPTED TO THE CRAPS ENVIRONMENT NOV/90/JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
      LOGICAL LDIAG
#include "mxcent.par"
#include "baslims.par"
      PARAMETER (D0 = 0.00 D00)
      DIMENSION WORK1(1)
      DIMENSION GRADKE(NCOORD),GRADFS(NCOORD),
     &          HESSKE(NCOORD,NCOORD),HESFS2(NCOORD,NCOORD)
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN
     &           *MXCONT*MXCONT)
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
C
      COMMON /GENCON/ NRCA,NRCB,CONTA(MXCONT),CONTB(MXCONT)
      LOGICAL SECOND, DOAX, DOAY, DOAZ, DOBX, DOBY, DOBZ
      COMMON/DF/DSHELL(MXAQNS), FSHELL(MXAQNS)
C
      DATA HALF, TWO /0.5D0,2.D0/
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTXOR(I,J) = XOR(I,J)
      ITRI(I,J) = MAX(I,J)*(MAX(I,J) - 1)/2 + MIN(I,J)
C
C     Get density and Fock elements
C
      IOFF=0
      ISKIP=NRCA*NRCB
      IF(LDIAG) ISKIP=NRCA*(NRCA+1)/2
      DO 200 IRCA=1,NRCA
       MAXB=NRCB
       IF(LDIAG) MAXB=IRCA
       DO 200 IRCB=1,MAXB
        IOFF=IOFF+1
        ICOMP=IOFF
        DO 100 IORBA = IDENA + IRCA, IDENA + NRCA*KHKTA,NRCA
        DO 150 IORBB = IDENB + IRCB, IDENB + NRCB*KHKTB,NRCB
         IORBAB = ITRI(IORBA,IORBB)
         DSHELL(ICOMP) = WORK1(KDEN-1+IORBAB)
         FSHELL(ICOMP) = WORK1(KFOC-1+IORBAB)
         ICOMP = ICOMP + ISKIP
150     CONTINUE
100    CONTINUE
       IF(LDIAG.AND.(IRCA.EQ.IRCB)) THEN
        CALL SSCAL(KHKTA*KHKTB,HALF,DSHELL(IOFF),ISKIP)
        CALL SSCAL(KHKTA,TWO,DSHELL(IOFF),ISKIP*(KHKTA+1))
        CALL SSCAL(KHKTA*KHKTB,HALF,FSHELL(IOFF),ISKIP)
        CALL SSCAL(KHKTA,TWO,FSHELL(IOFF),ISKIP*(KHKTA+1))
       ENDIF
200   CONTINUE
      MAXCMP = ISKIP*KHKTA*KHKTB
C
C     Multiply densities and Fock elements with integrals
C
c         IF (SECOND) THEN
c            DERTXX = DERTXX + DENS*TDERXX(ICOMP)
c            DERTXY = DERTXY + DENS*TDERXY(ICOMP)
c            DERTXZ = DERTXZ + DENS*TDERXZ(ICOMP)
c            DERTYY = DERTYY + DENS*TDERYY(ICOMP)
c            DERTYZ = DERTYZ + DENS*TDERYZ(ICOMP)
c            DERTZZ = DERTZZ + DENS*TDERZZ(ICOMP)
c            DERSXX = DERSXX - FOCK*SDERXX(ICOMP)
c            DERSXY = DERSXY - FOCK*SDERXY(ICOMP)
c            DERSXZ = DERSXZ - FOCK*SDERXZ(ICOMP)
c            DERSYY = DERSYY - FOCK*SDERYY(ICOMP)
c            DERSYZ = DERSYZ - FOCK*SDERYZ(ICOMP)
c            DERSZZ = DERSZZ - FOCK*SDERZZ(ICOMP)
c         END IF
c  200 CONTINUE
C
      DERT0 = SDOT(MAXCMP,DSHELL,1,WORK1(IT0000+1),1)
      DERTX = SDOT(MAXCMP,DSHELL,1,WORK1(IT000X+1),1)
      DERTY = SDOT(MAXCMP,DSHELL,1,WORK1(IT000Y+1),1)
      DERTZ = SDOT(MAXCMP,DSHELL,1,WORK1(IT000Z+1),1)
      DERSX = - SDOT(MAXCMP,FSHELL,1,WORK1(IS000X+1),1)
      DERSY = - SDOT(MAXCMP,FSHELL,1,WORK1(IS000Y+1),1)
      DERSZ = - SDOT(MAXCMP,FSHELL,1,WORK1(IS000Z+1),1)

CSSS      Write(6,*) "kinetic derivatives"
CSSS      Write(6,"(4(1x,f15.10))") (WORK1(IT000X+IOFF),ioff=1,Maxcmp)
CSSS      Write(6,"(4(1x,f15.10))") (WORK1(IT000Y+IOFF),ioff=1,Maxcmp)
CSSS      Write(6,"(4(1x,f15.10))") (WORK1(IT000Z+IOFF),ioff=1,Maxcmp)
CSSS      Write(6,*)
CSSS      Write(6,*) "Fock matrix element"
CSSS      Write(6,"(4(1x,f20.10))") (Fshell(IOFF), ioff=1,Maxcmp) 
CSSS      Write(6,*)
CSSS      Write(6,"(4(1x,f20.10))") DERTX, DERTY, DERTZ 
CSSS      Write(6,"(4(1x,f20.10))") DERSX, DERSY, DERSZ 
CSSS      Write(6,*) 

      IF(SECOND) THEN
       DERTXX = SDOT(MAXCMP,DSHELL,1,WORK1(IT00XX+1),1)
       DERTXY = SDOT(MAXCMP,DSHELL,1,WORK1(IT00XY+1),1)
       DERTXZ = SDOT(MAXCMP,DSHELL,1,WORK1(IT00XZ+1),1)
       DERTYY = SDOT(MAXCMP,DSHELL,1,WORK1(IT00YY+1),1)
       DERTYZ = SDOT(MAXCMP,DSHELL,1,WORK1(IT00YZ+1),1)
       DERTZZ = SDOT(MAXCMP,DSHELL,1,WORK1(IT00ZZ+1),1)
       DERSXX = - SDOT(MAXCMP,FSHELL,1,WORK1(IS00XX+1),1)
       DERSXY = - SDOT(MAXCMP,FSHELL,1,WORK1(IS00XY+1),1)
       DERSXZ = - SDOT(MAXCMP,FSHELL,1,WORK1(IS00XZ+1),1)
       DERSYY = - SDOT(MAXCMP,FSHELL,1,WORK1(IS00YY+1),1)
       DERSYZ = - SDOT(MAXCMP,FSHELL,1,WORK1(IS00YZ+1),1)
       DERSZZ = - SDOT(MAXCMP,FSHELL,1,WORK1(IS00ZZ+1),1)
      ENDIF

C
      NAX = 3*NCENTA - 2
      NAY = 3*NCENTA - 1
      NAZ = 3*NCENTA
      NBX = 3*NCENTB - 2
      NBY = 3*NCENTB - 1
      NBZ = 3*NCENTB
C
C     Undifferentiated kinetic energy
C     Note: One-center terms are not included!
C
      ENERKE = ENERKE + DERT0
C
C     ********************
C     ***** Gradient *****
C     ********************
C
      IAX  = IPTCNT(NAX,0)
      IAY  = IPTCNT(NAY,0)
      IAZ  = IPTCNT(NAZ,0)
      IBX  = IPTCNT(NBX,0)
      IBY  = IPTCNT(NBY,0)
      IBZ  = IPTCNT(NBZ,0)
      DOAX = IAX .NE. 0
      DOAY = IAY .NE. 0
      DOAZ = IAZ .NE. 0
      DOBX = IBX .NE. 0
      DOBY = IBY .NE. 0
      DOBZ = IBZ .NE. 0
C
C     A gradient elements:
C
      IF (DOAX) THEN
         GRADKE(IAX)   = GRADKE(IAX)   + DERTX
         GRADFS(IAX)   = GRADFS(IAX)   + DERSX
      END IF
      IF (DOAY) THEN
         GRADKE(IAY)   = GRADKE(IAY)   + DERTY
         GRADFS(IAY)   = GRADFS(IAY)   + DERSY
      END IF
      IF (DOAZ) THEN
         GRADKE(IAZ)   = GRADKE(IAZ)   + DERTZ
         GRADFS(IAZ)   = GRADFS(IAZ)   + DERSZ
CSSS         Write(6,"((1x,f20.10))") GRADKE(IAZ)
CSSS         Write(6,"((1x,f20.10))") GRADFS(IAZ)
      END IF
CSSS      Write(6,*) 
C
C     B gradient elements:
C
      IF (DOBX) THEN
         GRADKE(IBX)   = GRADKE(IBX)   - SBX*DERTX
         GRADFS(IBX)   = GRADFS(IBX)   - SBX*DERSX
      END IF
      IF (DOBY) THEN
         GRADKE(IBY)   = GRADKE(IBY)   - SBY*DERTY
         GRADFS(IBY)   = GRADFS(IBY)   - SBY*DERSY
      END IF
      IF (DOBZ) THEN
         GRADKE(IBZ)   = GRADKE(IBZ)   - SBZ*DERTZ
         GRADFS(IBZ)   = GRADFS(IBZ)   - SBZ*DERSZ
CSSS         Write(6,"((1x,f20.10))") GRADKE(IBZ)
CSSS         Write(6,"((1x,f20.10))") GRADFS(IBZ)
      END IF
CSSS      Write(6,*)
C
C     *******************
C     ***** Hessian *****
C     *******************
C
      IF (SECOND) THEN
       DO 300 IREP = 0, MAXLOP
        CHI = PT(IBTAND(ISYMOP,IREP))
        CSBX = CHI*SBX
        CSBY = CHI*SBY
        CSBZ = CHI*SBZ
        IF (NCENTA .EQ. NCENTB) THEN
         CASBX = CSBX + CSBX
         CASBY = CSBY + CSBY
         CASBZ = CSBZ + CSBZ
        ELSE
         CASBX = CSBX
         CASBY = CSBY
         CASBZ = CSBZ
        END IF
        IAX  = IPTCNT(NAX,IREP)
        IAY  = IPTCNT(NAY,IREP)
        IAZ  = IPTCNT(NAZ,IREP)
        IBX  = IPTCNT(NBX,IREP)
        IBY  = IPTCNT(NBY,IREP)
        IBZ  = IPTCNT(NBZ,IREP)
        DOAX = IAX .NE. 0
        DOAY = IAY .NE. 0
        DOAZ = IAZ .NE. 0
        DOBX = IBX .NE. 0
        DOBY = IBY .NE. 0
        DOBZ = IBZ .NE. 0
C
C           A-A and A-B Hessian elements:
C
C           First is Ax
C
        IF (DOAX) THEN
         HESSKE(IAX,IAX) = HESSKE(IAX,IAX) + DERTXX
         HESFS2(IAX,IAX) = HESFS2(IAX,IAX) + DERSXX
         IF (DOAY) THEN
          HESSKE(IAX,IAY) = HESSKE(IAX,IAY) + DERTXY
          HESFS2(IAX,IAY) = HESFS2(IAX,IAY) + DERSXY
         END IF
         IF (DOAZ) THEN
          HESSKE(IAX,IAZ) = HESSKE(IAX,IAZ) + DERTXZ
          HESFS2(IAX,IAZ) = HESFS2(IAX,IAZ) + DERSXZ
         END IF
         IF (DOBX) THEN
          HESSKE(IAX,IBX) = HESSKE(IAX,IBX) - CASBX*DERTXX
          HESFS2(IAX,IBX) = HESFS2(IAX,IBX) - CASBX*DERSXX
         END IF
         IF (DOBY) THEN
          HESSKE(IAX,IBY) = HESSKE(IAX,IBY) - CSBY*DERTXY
          HESFS2(IAX,IBY) = HESFS2(IAX,IBY) - CSBY*DERSXY
         END IF
         IF (DOBZ) THEN
          HESSKE(IAX,IBZ) = HESSKE(IAX,IBZ) - CSBZ*DERTXZ
          HESFS2(IAX,IBZ) = HESFS2(IAX,IBZ) - CSBZ*DERSXZ
         END IF
        END IF
C
C           First is Ay
C
        IF (DOAY) THEN
         HESSKE(IAY,IAY) = HESSKE(IAY,IAY) + DERTYY
         HESFS2(IAY,IAY) = HESFS2(IAY,IAY) + DERSYY
         IF (DOAZ) THEN
          HESSKE(IAY,IAZ) = HESSKE(IAY,IAZ) + DERTYZ
          HESFS2(IAY,IAZ) = HESFS2(IAY,IAZ) + DERSYZ
         END IF
         IF (DOBX) THEN
          HESSKE(IAY,IBX) = HESSKE(IAY,IBX) - CSBX*DERTXY
          HESFS2(IAY,IBX) = HESFS2(IAY,IBX) - CSBX*DERSXY
         END IF
         IF (DOBY) THEN
          HESSKE(IAY,IBY) = HESSKE(IAY,IBY) - CASBY*DERTYY
          HESFS2(IAY,IBY) = HESFS2(IAY,IBY) - CASBY*DERSYY
         END IF
         IF (DOBZ) THEN
          HESSKE(IAY,IBZ) = HESSKE(IAY,IBZ) - CSBZ*DERTYZ
          HESFS2(IAY,IBZ) = HESFS2(IAY,IBZ) - CSBZ*DERSYZ
         END IF
        END IF
C
C           First is Az
C
        IF (DOAZ) THEN
         HESSKE(IAZ,IAZ) = HESSKE(IAZ,IAZ) + DERTZZ
         HESFS2(IAZ,IAZ) = HESFS2(IAZ,IAZ) + DERSZZ
         IF (DOBX) THEN
          HESSKE(IAZ,IBX) = HESSKE(IAZ,IBX) - CSBX*DERTXZ
          HESFS2(IAZ,IBX) = HESFS2(IAZ,IBX) - CSBX*DERSXZ
         END IF
         IF (DOBY) THEN
          HESSKE(IAZ,IBY) = HESSKE(IAZ,IBY) - CSBY*DERTYZ
          HESFS2(IAZ,IBY) = HESFS2(IAZ,IBY) - CSBY*DERSYZ
         END IF
         IF (DOBZ) THEN
          HESSKE(IAZ,IBZ) = HESSKE(IAZ,IBZ) - CASBZ*DERTZZ
          HESFS2(IAZ,IBZ) = HESFS2(IAZ,IBZ) - CASBZ*DERSZZ
         END IF
        END IF
C
C           B-B Hessian elements:
C
C           First is Bx
C
        IF (DOBX) THEN
         HESSKE(IBX,IBX) = HESSKE(IBX,IBX) + DERTXX
         HESFS2(IBX,IBX) = HESFS2(IBX,IBX) + DERSXX
         IF (DOBY) THEN
          HESSKE(IBX,IBY) = HESSKE(IBX,IBY) + CSBX*CSBY*DERTXY
          HESFS2(IBX,IBY) = HESFS2(IBX,IBY) + CSBX*CSBY*DERSXY
         END IF
         IF (DOBZ) THEN
          HESSKE(IBX,IBZ) = HESSKE(IBX,IBZ) + CSBX*CSBZ*DERTXZ
          HESFS2(IBX,IBZ) = HESFS2(IBX,IBZ) + CSBX*CSBZ*DERSXZ
         END IF
        END IF
C
C           First is By
C
        IF (DOBY) THEN
         HESSKE(IBY,IBY) = HESSKE(IBY,IBY) + DERTYY
         HESFS2(IBY,IBY) = HESFS2(IBY,IBY) + DERSYY
         IF (DOBZ) THEN
          HESSKE(IBY,IBZ) = HESSKE(IBY,IBZ) + CSBY*CSBZ*DERTYZ
          HESFS2(IBY,IBZ) = HESFS2(IBY,IBZ) + CSBY*CSBZ*DERSYZ
         END IF
        END IF
C
C           First is Bz
C
        IF (DOBZ) THEN
         HESSKE(IBZ,IBZ) = HESSKE(IBZ,IBZ) + DERTZZ
         HESFS2(IBZ,IBZ) = HESFS2(IBZ,IBZ) + DERSZZ
        END IF
300    CONTINUE
      END IF
      RETURN
      END
