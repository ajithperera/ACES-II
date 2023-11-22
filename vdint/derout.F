      SUBROUTINE DEROUT(NINTYP,NCENT1,NCENT2,NCENT3,NCENT4,
     &                  ISOA,ISOB,ISOC,ISOD,
     &                  SAX,SAY,SAZ,SBX,SBY,SBZ,
     &                  SCX,SCY,SCZ,SDX,SDY,SDZ,
     &                  DERONE,DERTWO,TWOCEN,THRCEN,FOUCEN,
     &                  IPRINT,DERIV,BIGVEC,ENEREE,GRADEE,
     &                  HESSEE,NCOORD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "mxcent.par"
#include "baslims.par"
      PARAMETER (D1 = 1.00D00, D2 = 2.00D00)
      LOGICAL TWOCEN, THRCEN, DERONE, DERTWO, BIGVEC,
     *        DOAX, DOAY, DOAZ, DOBX, DOBY, DOBZ,
     *        DOCX, DOCY, DOCZ, DODX, DODY, DODZ, FOUCEN
      DIMENSION DERIV(54)
C
      DIMENSION GRADEE(NCOORD),HESSEE(NCOORD,NCOORD)
      CHARACTER NAMEX*6
      LOGICAL DCORD, DCORGD, NOORBT, DOPERT
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
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
      IF (IPRINT .GE.10) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 1110) TWOCEN, THRCEN
         WRITE (LUPRI, 1120) NCENT1, NCENT2, NCENT3, NCENT4
      END IF
C
C     ***************************
C     ***** Two-Center Case *****
C     ***************************
C
      IF (TWOCEN)THEN
         NAX = 3*NCENT1 - 2
         NAY = 3*NCENT1 - 1
         NAZ = 3*NCENT1
         NBX = 3*NCENT2 - 2
         NBY = 3*NCENT2 - 1
         NBZ = 3*NCENT2
C
C        ***** Gradient *****
C
         IF (DERONE) THEN
            DAX = DERIV(1)
            DAY = DERIV(2)
            DAZ = DERIV(3)
            IAX = IPTCNT(NAX,0)
            IAY = IPTCNT(NAY,0)
            IAZ = IPTCNT(NAZ,0)
            IBX = IPTCNT(NBX,0)
            IBY = IPTCNT(NBY,0)
            IBZ = IPTCNT(NBZ,0)
            DOAX = IAX .NE. 0
            DOAY = IAY .NE. 0
            DOAZ = IAZ .NE. 0
            DOBX = IBX .NE. 0
            DOBY = IBY .NE. 0
            DOBZ = IBZ .NE. 0
C
C           A electron-repulsion gradient elements:
C
            IF (DOAX) GRADEE(IAX) = GRADEE(IAX) + SAX*DAX
            IF (DOAY) GRADEE(IAY) = GRADEE(IAY) + SAY*DAY
            IF (DOAZ) GRADEE(IAZ) = GRADEE(IAZ) + SAZ*DAZ
C
C           B electron-repulsion gradient elements:
C
            IF (DOBX) GRADEE(IBX) = GRADEE(IBX) - SBX*DAX
            IF (DOBY) GRADEE(IBY) = GRADEE(IBY) - SBY*DAY
            IF (DOBZ) GRADEE(IBZ) = GRADEE(IBZ) - SBZ*DAZ
         END IF
C
C        ***** Hessian *****
C
         IF (DERTWO) THEN
            FAB = D1
            IF (NCENT1 .EQ. NCENT2) FAB = D2
C
            DIAAXX = DERIV(4)
            DIAAXY = DERIV(5)
            DIAAXZ = DERIV(6)
            DIAAYY = DERIV(7)
            DIAAYZ = DERIV(8)
            DIAAZZ = DERIV(9)
C
            DAAXX =         DIAAXX
            DAAXY = SAX*SAY*DIAAXY
            DAAXZ = SAX*SAZ*DIAAXZ
            DAAYY =         DIAAYY
            DAAYZ = SAY*SAZ*DIAAYZ
            DAAZZ =         DIAAZZ
C
            DABXX = - FAB*SAX*SBX*DIAAXX
            DABXY = -     SAX*SBY*DIAAXY
            DABXZ = -     SAX*SBZ*DIAAXZ
            DABYX = -     SAY*SBX*DIAAXY
            DABYY = - FAB*SAY*SBY*DIAAYY
            DABYZ = -     SAY*SBZ*DIAAYZ
            DABZX = -     SAZ*SBX*DIAAXZ
            DABZY = -     SAZ*SBY*DIAAYZ
            DABZZ = - FAB*SAZ*SBZ*DIAAZZ
C
            DBBXX =         DIAAXX
            DBBXY = SBX*SBY*DIAAXY
            DBBXZ = SBX*SBZ*DIAAXZ
            DBBYY =         DIAAYY
            DBBYZ = SBY*SBZ*DIAAYZ
            DBBZZ =         DIAAZZ
C
C           Loop over irreps
C
            DO 100 IREP = 0, MAXLOP
               CHIAB = PT(IBTAND(ISOA,IREP))*PT(IBTAND(ISOB,IREP))
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
C              A-A electron-repulsion Hessian elements:
C
               IF (DOAX)
     *            HESSEE(IAX,IAX) = HESSEE(IAX,IAX) + DAAXX
               IF (DOAX.AND.DOAY)
     *            HESSEE(IAX,IAY) = HESSEE(IAX,IAY) + DAAXY
               IF (DOAX.AND.DOAZ)
     *            HESSEE(IAX,IAZ) = HESSEE(IAX,IAZ) + DAAXZ
               IF (DOAY)
     *            HESSEE(IAY,IAY) = HESSEE(IAY,IAY) + DAAYY
               IF (DOAY.AND.DOAZ)
     *            HESSEE(IAY,IAZ) = HESSEE(IAY,IAZ) + DAAYZ
               IF (DOAZ)
     *            HESSEE(IAZ,IAZ) = HESSEE(IAZ,IAZ) + DAAZZ
C
C              A-B electron-repulsion Hessian elements:
C
               IF (DOAX.AND.DOBX)
     *            HESSEE(IAX,IBX) = HESSEE(IAX,IBX) + CHIAB*DABXX
               IF (DOAX.AND.DOBY)
     *            HESSEE(IAX,IBY) = HESSEE(IAX,IBY) + CHIAB*DABXY
               IF (DOAX.AND.DOBZ)
     *            HESSEE(IAX,IBZ) = HESSEE(IAX,IBZ) + CHIAB*DABXZ
               IF (DOAY.AND.DOBX)
     *            HESSEE(IAY,IBX) = HESSEE(IAY,IBX) + CHIAB*DABYX
               IF (DOAY.AND.DOBY)
     *            HESSEE(IAY,IBY) = HESSEE(IAY,IBY) + CHIAB*DABYY
               IF (DOAY.AND.DOBZ)
     *            HESSEE(IAY,IBZ) = HESSEE(IAY,IBZ) + CHIAB*DABYZ
               IF (DOAZ.AND.DOBX)
     *            HESSEE(IAZ,IBX) = HESSEE(IAZ,IBX) + CHIAB*DABZX
               IF (DOAZ.AND.DOBY)
     *            HESSEE(IAZ,IBY) = HESSEE(IAZ,IBY) + CHIAB*DABZY
               IF (DOAZ.AND.DOBZ)
     *            HESSEE(IAZ,IBZ) = HESSEE(IAZ,IBZ) + CHIAB*DABZZ
C
C              B-B electron-repulsion Hessian elements:
C
               IF (DOBX)
     *            HESSEE(IBX,IBX) = HESSEE(IBX,IBX) + DBBXX
               IF (DOBX.AND.DOBY)
     *            HESSEE(IBX,IBY) = HESSEE(IBX,IBY) + DBBXY
               IF (DOBX.AND.DOBZ)
     *            HESSEE(IBX,IBZ) = HESSEE(IBX,IBZ) + DBBXZ
               IF (DOBY)
     *            HESSEE(IBY,IBY) = HESSEE(IBY,IBY) + DBBYY
               IF (DOBY.AND.DOBZ)
     *            HESSEE(IBY,IBZ) = HESSEE(IBY,IBZ) + DBBYZ
               IF (DOBZ)
     *            HESSEE(IBZ,IBZ) = HESSEE(IBZ,IBZ) + DBBZZ
  100       CONTINUE
         END IF
C
C     *****************************
C     ***** Three-Center Case *****
C     *****************************
C
      ELSE IF (THRCEN) THEN
         NAX = 3*NCENT1 - 2
         NAY = 3*NCENT1 - 1
         NAZ = 3*NCENT1
         NBX = 3*NCENT2 - 2
         NBY = 3*NCENT2 - 1
         NBZ = 3*NCENT2
         NCX = 3*NCENT3 - 2
         NCY = 3*NCENT3 - 1
         NCZ = 3*NCENT3
C
C        ***** Gradient *****
C
         IF (DERONE) THEN
            DAX = DERIV(1)
            DAY = DERIV(2)
            DAZ = DERIV(3)
            DBX = DERIV(4)
            DBY = DERIV(5)
            DBZ = DERIV(6)
            IAX = IPTCNT(NAX,0)
            IAY = IPTCNT(NAY,0)
            IAZ = IPTCNT(NAZ,0)
            IBX = IPTCNT(NBX,0)
            IBY = IPTCNT(NBY,0)
            IBZ = IPTCNT(NBZ,0)
            ICX = IPTCNT(NCX,0)
            ICY = IPTCNT(NCY,0)
            ICZ = IPTCNT(NCZ,0)
            DOAX = IAX .NE. 0
            DOAY = IAY .NE. 0
            DOAZ = IAZ .NE. 0
            DOBX = IBX .NE. 0
            DOBY = IBY .NE. 0
            DOBZ = IBZ .NE. 0
            DOCX = ICX .NE. 0
            DOCY = ICY .NE. 0
            DOCZ = ICZ .NE. 0
C
C           A electron-repulsion gradient elements:
C
            IF (DOAX) GRADEE(IAX) = GRADEE(IAX) + SAX*DAX
            IF (DOAY) GRADEE(IAY) = GRADEE(IAY) + SAY*DAY
            IF (DOAZ) GRADEE(IAZ) = GRADEE(IAZ) + SAZ*DAZ
C
C           B electron-repulsion gradient elements:
C
            IF (DOBX) GRADEE(IBX) = GRADEE(IBX) + SBX*DBX
            IF (DOBY) GRADEE(IBY) = GRADEE(IBY) + SBY*DBY
            IF (DOBZ) GRADEE(IBZ) = GRADEE(IBZ) + SBZ*DBZ
C
C           C electron-repulsion gradient elements:
C
            IF (DOCX) GRADEE(ICX) = GRADEE(ICX) - SCX*(DAX + DBX)
            IF (DOCY) GRADEE(ICY) = GRADEE(ICY) - SCY*(DAY + DBY)
            IF (DOCZ) GRADEE(ICZ) = GRADEE(ICZ) - SCZ*(DAZ + DBZ)
         END IF
C
C        ***** Hessian *****
C
         IF (DERTWO) THEN
            FAB = D1
            FAC = D1
            FBC = D1
            IF (NCENT1 .EQ. NCENT2) FAB = D2
            IF (NCENT1 .EQ. NCENT3) FAC = D2
            IF (NCENT2 .EQ. NCENT3) FBC = D2
C
            DIAAXX = DERIV(7)
            DIAAXY = DERIV(8)
            DIAAXZ = DERIV(9)
            DIAAYY = DERIV(10)
            DIAAYZ = DERIV(11)
            DIAAZZ = DERIV(12)
C
            DIBBXX = DERIV(13)
            DIBBXY = DERIV(14)
            DIBBXZ = DERIV(15)
            DIBBYY = DERIV(16)
            DIBBYZ = DERIV(17)
            DIBBZZ = DERIV(18)
C
            DIABXX = DERIV(19)
            DIABXY = DERIV(20)
            DIABXZ = DERIV(21)
            DIABYX = DERIV(22)
            DIABYY = DERIV(23)
            DIABYZ = DERIV(24)
            DIABZX = DERIV(25)
            DIABZY = DERIV(26)
            DIABZZ = DERIV(27)
C
            DAAXX =         DIAAXX
            DAAXY = SAX*SAY*DIAAXY
            DAAXZ = SAX*SAZ*DIAAXZ
            DAAYY =         DIAAYY
            DAAYZ = SAY*SAZ*DIAAYZ
            DAAZZ =         DIAAZZ
C
            DBBXX =         DIBBXX
            DBBXY = SBX*SBY*DIBBXY
            DBBXZ = SBX*SBZ*DIBBXZ
            DBBYY =         DIBBYY
            DBBYZ = SBY*SBZ*DIBBYZ
            DBBZZ =         DIBBZZ
C
            DABXX = FAB*SAX*SBX*DIABXX
            DABXY =     SAX*SBY*DIABXY
            DABXZ =     SAX*SBZ*DIABXZ
            DABYX =     SAY*SBX*DIABYX
            DABYY = FAB*SAY*SBY*DIABYY
            DABYZ =     SAY*SBZ*DIABYZ
            DABZX =     SAZ*SBX*DIABZX
            DABZY =     SAZ*SBY*DIABZY
            DABZZ = FAB*SAZ*SBZ*DIABZZ
C
            DACXX = - FAC*SAX*SCX*(DIAAXX + DIABXX)
            DACXY = -     SAX*SCY*(DIAAXY + DIABXY)
            DACXZ = -     SAX*SCZ*(DIAAXZ + DIABXZ)
            DACYX = -     SAY*SCX*(DIAAXY + DIABYX)
            DACYY = - FAC*SAY*SCY*(DIAAYY + DIABYY)
            DACYZ = -     SAY*SCZ*(DIAAYZ + DIABYZ)
            DACZX = -     SAZ*SCX*(DIAAXZ + DIABZX)
            DACZY = -     SAZ*SCY*(DIAAYZ + DIABZY)
            DACZZ = - FAC*SAZ*SCZ*(DIAAZZ + DIABZZ)
C
            DBCXX = - FBC*SBX*SCX*(DIBBXX + DIABXX)
            DBCXY = -     SBX*SCY*(DIBBXY + DIABYX)
            DBCXZ = -     SBX*SCZ*(DIBBXZ + DIABZX)
            DBCYX = -     SBY*SCX*(DIBBXY + DIABXY)
            DBCYY = - FBC*SBY*SCY*(DIBBYY + DIABYY)
            DBCYZ = -     SBY*SCZ*(DIBBYZ + DIABZY)
            DBCZX = -     SBZ*SCX*(DIBBXZ + DIABXZ)
            DBCZY = -     SBZ*SCY*(DIBBYZ + DIABYZ)
            DBCZZ = - FBC*SBZ*SCZ*(DIBBZZ + DIABZZ)
C
            DCCXX =           (DIAAXX + DIABXX + DIABXX + DIBBXX)
            DCCXY =   SCX*SCY*(DIAAXY + DIABXY + DIABYX + DIBBXY)
            DCCXZ =   SCX*SCZ*(DIAAXZ + DIABXZ + DIABZX + DIBBXZ)
            DCCYY =           (DIAAYY + DIABYY + DIABYY + DIBBYY)
            DCCYZ =   SCY*SCZ*(DIAAYZ + DIABYZ + DIABZY + DIBBYZ)
            DCCZZ =           (DIAAZZ + DIABZZ + DIABZZ + DIBBZZ)
C
C           Loop over irreps
C
            DO 200 IREP = 0, MAXLOP
               CHIA = PT(IBTAND(ISOA,IREP))
               CHIB = PT(IBTAND(ISOB,IREP))
               CHIC = PT(IBTAND(ISOC,IREP))
               CHIAB = CHIA*CHIB
               CHIAC = CHIA*CHIC
               CHIBC = CHIB*CHIC
               IAX  = IPTCNT(NAX,IREP)
               IAY  = IPTCNT(NAY,IREP)
               IAZ  = IPTCNT(NAZ,IREP)
               IBX  = IPTCNT(NBX,IREP)
               IBY  = IPTCNT(NBY,IREP)
               IBZ  = IPTCNT(NBZ,IREP)
               ICX  = IPTCNT(NCX,IREP)
               ICY  = IPTCNT(NCY,IREP)
               ICZ  = IPTCNT(NCZ,IREP)
               DOAX = IAX .NE. 0
               DOAY = IAY .NE. 0
               DOAZ = IAZ .NE. 0
               DOBX = IBX .NE. 0
               DOBY = IBY .NE. 0
               DOBZ = IBZ .NE. 0
               DOCX = ICX .NE. 0
               DOCY = ICY .NE. 0
               DOCZ = ICZ .NE. 0
C
C              A-A electron-repulsion Hessian elements:
C
               IF (DOAX)
     *            HESSEE(IAX,IAX) = HESSEE(IAX,IAX) + DAAXX
               IF (DOAX.AND.DOAY)
     *            HESSEE(IAX,IAY) = HESSEE(IAX,IAY) + DAAXY
               IF (DOAX.AND.DOAZ)
     *            HESSEE(IAX,IAZ) = HESSEE(IAX,IAZ) + DAAXZ
               IF (DOAY)
     *            HESSEE(IAY,IAY) = HESSEE(IAY,IAY) + DAAYY
               IF (DOAY.AND.DOAZ)
     *            HESSEE(IAY,IAZ) = HESSEE(IAY,IAZ) + DAAYZ
               IF (DOAZ)
     *            HESSEE(IAZ,IAZ) = HESSEE(IAZ,IAZ) + DAAZZ
C
C              A-B electron-repulsion Hessian elements:
C
               IF (DOAX.AND.DOBX)
     *            HESSEE(IAX,IBX) = HESSEE(IAX,IBX) + CHIAB*DABXX
               IF (DOAX.AND.DOBY)
     *            HESSEE(IAX,IBY) = HESSEE(IAX,IBY) + CHIAB*DABXY
               IF (DOAX.AND.DOBZ)
     *            HESSEE(IAX,IBZ) = HESSEE(IAX,IBZ) + CHIAB*DABXZ
               IF (DOAY.AND.DOBX)
     *            HESSEE(IAY,IBX) = HESSEE(IAY,IBX) + CHIAB*DABYX
               IF (DOAY.AND.DOBY)
     *            HESSEE(IAY,IBY) = HESSEE(IAY,IBY) + CHIAB*DABYY
               IF (DOAY.AND.DOBZ)
     *            HESSEE(IAY,IBZ) = HESSEE(IAY,IBZ) + CHIAB*DABYZ
               IF (DOAZ.AND.DOBX)
     *            HESSEE(IAZ,IBX) = HESSEE(IAZ,IBX) + CHIAB*DABZX
               IF (DOAZ.AND.DOBY)
     *            HESSEE(IAZ,IBY) = HESSEE(IAZ,IBY) + CHIAB*DABZY
               IF (DOAZ.AND.DOBZ)
     *            HESSEE(IAZ,IBZ) = HESSEE(IAZ,IBZ) + CHIAB*DABZZ
C
C              A-C electron-repulsion Hessian elements:
C
               IF (DOAX.AND.DOCX)
     *            HESSEE(IAX,ICX) = HESSEE(IAX,ICX) + CHIAC*DACXX
               IF (DOAX.AND.DOCY)
     *            HESSEE(IAX,ICY) = HESSEE(IAX,ICY) + CHIAC*DACXY
               IF (DOAX.AND.DOCZ)
     *            HESSEE(IAX,ICZ) = HESSEE(IAX,ICZ) + CHIAC*DACXZ
               IF (DOAY.AND.DOCX)
     *            HESSEE(IAY,ICX) = HESSEE(IAY,ICX) + CHIAC*DACYX
               IF (DOAY.AND.DOCY)
     *            HESSEE(IAY,ICY) = HESSEE(IAY,ICY) + CHIAC*DACYY
               IF (DOAY.AND.DOCZ)
     *            HESSEE(IAY,ICZ) = HESSEE(IAY,ICZ) + CHIAC*DACYZ
               IF (DOAZ.AND.DOCX)
     *            HESSEE(IAZ,ICX) = HESSEE(IAZ,ICX) + CHIAC*DACZX
               IF (DOAZ.AND.DOCY)
     *            HESSEE(IAZ,ICY) = HESSEE(IAZ,ICY) + CHIAC*DACZY
               IF (DOAZ.AND.DOCZ)
     *            HESSEE(IAZ,ICZ) = HESSEE(IAZ,ICZ) + CHIAC*DACZZ
C
C              B-B electron-repulsion Hessian elements:
C
               IF (DOBX)
     *            HESSEE(IBX,IBX) = HESSEE(IBX,IBX) + DBBXX
               IF (DOBX.AND.DOBY)
     *            HESSEE(IBX,IBY) = HESSEE(IBX,IBY) + DBBXY
               IF (DOBX.AND.DOBZ)
     *            HESSEE(IBX,IBZ) = HESSEE(IBX,IBZ) + DBBXZ
               IF (DOBY)
     *            HESSEE(IBY,IBY) = HESSEE(IBY,IBY) + DBBYY
               IF (DOBY.AND.DOBZ)
     *            HESSEE(IBY,IBZ) = HESSEE(IBY,IBZ) + DBBYZ
               IF (DOBZ)
     *            HESSEE(IBZ,IBZ) = HESSEE(IBZ,IBZ) + DBBZZ
C
C              B-C electron-repulsion Hessian elements:
C
               IF (DOBX.AND.DOCX)
     *            HESSEE(IBX,ICX) = HESSEE(IBX,ICX) + CHIBC*DBCXX
               IF (DOBX.AND.DOCY)
     *            HESSEE(IBX,ICY) = HESSEE(IBX,ICY) + CHIBC*DBCXY
               IF (DOBX.AND.DOCZ)
     *            HESSEE(IBX,ICZ) = HESSEE(IBX,ICZ) + CHIBC*DBCXZ
               IF (DOBY.AND.DOCX)
     *            HESSEE(IBY,ICX) = HESSEE(IBY,ICX) + CHIBC*DBCYX
               IF (DOBY.AND.DOCY)
     *            HESSEE(IBY,ICY) = HESSEE(IBY,ICY) + CHIBC*DBCYY
               IF (DOBY.AND.DOCZ)
     *            HESSEE(IBY,ICZ) = HESSEE(IBY,ICZ) + CHIBC*DBCYZ
               IF (DOBZ.AND.DOCX)
     *            HESSEE(IBZ,ICX) = HESSEE(IBZ,ICX) + CHIBC*DBCZX
               IF (DOBZ.AND.DOCY)
     *            HESSEE(IBZ,ICY) = HESSEE(IBZ,ICY) + CHIBC*DBCZY
               IF (DOBZ.AND.DOCZ)
     *            HESSEE(IBZ,ICZ) = HESSEE(IBZ,ICZ) + CHIBC*DBCZZ
C
C              C-C electron-repulsion Hessian elements:
C
               IF (DOCX)
     *            HESSEE(ICX,ICX) = HESSEE(ICX,ICX) + DCCXX
               IF (DOCX.AND.DOCY)
     *            HESSEE(ICX,ICY) = HESSEE(ICX,ICY) + DCCXY
               IF (DOCX.AND.DOCZ)
     *            HESSEE(ICX,ICZ) = HESSEE(ICX,ICZ) + DCCXZ
               IF (DOCY)
     *            HESSEE(ICY,ICY) = HESSEE(ICY,ICY) + DCCYY
               IF (DOCY.AND.DOCZ)
     *            HESSEE(ICY,ICZ) = HESSEE(ICY,ICZ) + DCCYZ
               IF (DOCZ)
     *            HESSEE(ICZ,ICZ) = HESSEE(ICZ,ICZ) + DCCZZ
  200       CONTINUE
         END IF
C
C     ****************************
C     ***** Four-Center Case *****
C     ****************************
C
       ELSE
         NAX = 3*NCENT1 - 2
         NAY = 3*NCENT1 - 1
         NAZ = 3*NCENT1
         NBX = 3*NCENT2 - 2
         NBY = 3*NCENT2 - 1
         NBZ = 3*NCENT2
         NCX = 3*NCENT3 - 2
         NCY = 3*NCENT3 - 1
         NCZ = 3*NCENT3
         NDX = 3*NCENT4 - 2
         NDY = 3*NCENT4 - 1
         NDZ = 3*NCENT4
C
C        ***** Gradient *****
C
         IF (DERONE) THEN
            DAX = DERIV(1)
            DAY = DERIV(2)
            DAZ = DERIV(3)
            DBX = DERIV(4)
            DBY = DERIV(5)
            DBZ = DERIV(6)
            DCX = DERIV(7)
            DCY = DERIV(8)
            DCZ = DERIV(9)
            IAX = IPTCNT(NAX,0)
            IAY = IPTCNT(NAY,0)
            IAZ = IPTCNT(NAZ,0)
            IBX = IPTCNT(NBX,0)
            IBY = IPTCNT(NBY,0)
            IBZ = IPTCNT(NBZ,0)
            ICX = IPTCNT(NCX,0)
            ICY = IPTCNT(NCY,0)
            ICZ = IPTCNT(NCZ,0)
            IDX = IPTCNT(NDX,0)
            IDY = IPTCNT(NDY,0)
            IDZ = IPTCNT(NDZ,0)
            DOAX = IAX .NE. 0
            DOAY = IAY .NE. 0
            DOAZ = IAZ .NE. 0
            DOBX = IBX .NE. 0
            DOBY = IBY .NE. 0
            DOBZ = IBZ .NE. 0
            DOCX = ICX .NE. 0
            DOCY = ICY .NE. 0
            DOCZ = ICZ .NE. 0
            DODX = IDX .NE. 0
            DODY = IDY .NE. 0
            DODZ = IDZ .NE. 0
C
C           A electron-repulsion gradient elements:
C
            IF (DOAX) GRADEE(IAX) = GRADEE(IAX) + SAX*DAX
            IF (DOAY) GRADEE(IAY) = GRADEE(IAY) + SAY*DAY
            IF (DOAZ) GRADEE(IAZ) = GRADEE(IAZ) + SAZ*DAZ
C
C           B electron-repulsion gradient elements:
C
            IF (DOBX) GRADEE(IBX) = GRADEE(IBX) + SBX*DBX
            IF (DOBY) GRADEE(IBY) = GRADEE(IBY) + SBY*DBY
            IF (DOBZ) GRADEE(IBZ) = GRADEE(IBZ) + SBZ*DBZ
C
C           C electron-repulsion gradient elements:
C
            IF (DOCX) GRADEE(ICX) = GRADEE(ICX) + SCX*DCX
            IF (DOCY) GRADEE(ICY) = GRADEE(ICY) + SCY*DCY
            IF (DOCZ) GRADEE(ICZ) = GRADEE(ICZ) + SCZ*DCZ
C
C           D electron-repulsion gradient elements:
C
            IF (DODX) GRADEE(IDX) = GRADEE(IDX) - SDX*(DAX + DBX + DCX)
            IF (DODY) GRADEE(IDY) = GRADEE(IDY) - SDY*(DAY + DBY + DCY)
            IF (DODZ) GRADEE(IDZ) = GRADEE(IDZ) - SDZ*(DAZ + DBZ + DCZ)
         END IF
C
C        ***** Hessian *****
C
         IF (DERTWO) THEN
            FAB = D1
            FAC = D1
            FAD = D1
            FBC = D1
            FBD = D1
            FCD = D1
            IF (NCENT1 .EQ. NCENT2) FAB = D2
            IF (NCENT1 .EQ. NCENT3) FAC = D2
            IF (NCENT1 .EQ. NCENT4) FAD = D2
            IF (NCENT2 .EQ. NCENT3) FBC = D2
            IF (NCENT2 .EQ. NCENT4) FBD = D2
            IF (NCENT3 .EQ. NCENT4) FCD = D2
            IF (BIGVEC) THEN
               IF (IAX .EQ. IBX) FAB = D2
               IF (IAX .EQ. ICX) FAC = D2
               IF (IAX .EQ. IDX) FAD = D2
               IF (IBX .EQ. ICX) FBC = D2
               IF (IBX .EQ. IDX) FBD = D2
               IF (ICX .EQ. IDX) FCD = D2
               STOP 
            END IF
            DIAAXX = DERIV(10)
            DIAAXY = DERIV(11)
            DIAAXZ = DERIV(12)
            DIAAYY = DERIV(13)
            DIAAYZ = DERIV(14)
            DIAAZZ = DERIV(15)
C
            DIBBXX = DERIV(16)
            DIBBXY = DERIV(17)
            DIBBXZ = DERIV(18)
            DIBBYY = DERIV(19)
            DIBBYZ = DERIV(20)
            DIBBZZ = DERIV(21)
C
            DICCXX = DERIV(22)
            DICCXY = DERIV(23)
            DICCXZ = DERIV(24)
            DICCYY = DERIV(25)
            DICCYZ = DERIV(26)
            DICCZZ = DERIV(27)
C
            DIABXX = DERIV(28)
            DIABXY = DERIV(29)
            DIABXZ = DERIV(30)
            DIABYX = DERIV(31)
            DIABYY = DERIV(32)
            DIABYZ = DERIV(33)
            DIABZX = DERIV(34)
            DIABZY = DERIV(35)
            DIABZZ = DERIV(36)
C
            DIACXX = DERIV(37)
            DIACXY = DERIV(38)
            DIACXZ = DERIV(39)
            DIACYX = DERIV(40)
            DIACYY = DERIV(41)
            DIACYZ = DERIV(42)
            DIACZX = DERIV(43)
            DIACZY = DERIV(44)
            DIACZZ = DERIV(45)
C
            DIBCXX = DERIV(46)
            DIBCXY = DERIV(47)
            DIBCXZ = DERIV(48)
            DIBCYX = DERIV(49)
            DIBCYY = DERIV(50)
            DIBCYZ = DERIV(51)
            DIBCZX = DERIV(52)
            DIBCZY = DERIV(53)
            DIBCZZ = DERIV(54)
C
            DAAXX =         DIAAXX
            DAAXY = SAX*SAY*DIAAXY
            DAAXZ = SAX*SAZ*DIAAXZ
            DAAYY =         DIAAYY
            DAAYZ = SAY*SAZ*DIAAYZ
            DAAZZ =         DIAAZZ
C
            DBBXX =         DIBBXX
            DBBXY = SBX*SBY*DIBBXY
            DBBXZ = SBX*SBZ*DIBBXZ
            DBBYY =         DIBBYY
            DBBYZ = SBY*SBZ*DIBBYZ
            DBBZZ =         DIBBZZ
C
            DCCXX =         DICCXX
            DCCXY = SCX*SCY*DICCXY
            DCCXZ = SCX*SCZ*DICCXZ
            DCCYY =         DICCYY
            DCCYZ = SCY*SCZ*DICCYZ
            DCCZZ =         DICCZZ
C
            DABXX = FAB*SAX*SBX*DIABXX
            DABXY =     SAX*SBY*DIABXY
            DABXZ =     SAX*SBZ*DIABXZ
            DABYX =     SAY*SBX*DIABYX
            DABYY = FAB*SAY*SBY*DIABYY
            DABYZ =     SAY*SBZ*DIABYZ
            DABZX =     SAZ*SBX*DIABZX
            DABZY =     SAZ*SBY*DIABZY
            DABZZ = FAB*SAZ*SBZ*DIABZZ
C
            DACXX = FAC*SAX*SCX*DIACXX
            DACXY =     SAX*SCY*DIACXY
            DACXZ =     SAX*SCZ*DIACXZ
            DACYX =     SAY*SCX*DIACYX
            DACYY = FAC*SAY*SCY*DIACYY
            DACYZ =     SAY*SCZ*DIACYZ
            DACZX =     SAZ*SCX*DIACZX
            DACZY =     SAZ*SCY*DIACZY
            DACZZ = FAC*SAZ*SCZ*DIACZZ
C
            DBCXX = FBC*SBX*SCX*DIBCXX
            DBCXY =     SBX*SCY*DIBCXY
            DBCXZ =     SBX*SCZ*DIBCXZ
            DBCYX =     SBY*SCX*DIBCYX
            DBCYY = FBC*SBY*SCY*DIBCYY
            DBCYZ =     SBY*SCZ*DIBCYZ
            DBCZX =     SBZ*SCX*DIBCZX
            DBCZY =     SBZ*SCY*DIBCZY
            DBCZZ = FBC*SBZ*SCZ*DIBCZZ
C
            DADXX = - FAD*SAX*SDX*(DIAAXX + DIABXX + DIACXX)
            DADXY = -     SAX*SDY*(DIAAXY + DIABXY + DIACXY)
            DADXZ = -     SAX*SDZ*(DIAAXZ + DIABXZ + DIACXZ)
            DADYX = -     SAY*SDX*(DIAAXY + DIABYX + DIACYX)
            DADYY = - FAD*SAY*SDY*(DIAAYY + DIABYY + DIACYY)
            DADYZ = -     SAY*SDZ*(DIAAYZ + DIABYZ + DIACYZ)
            DADZX = -     SAZ*SDX*(DIAAXZ + DIABZX + DIACZX)
            DADZY = -     SAZ*SDY*(DIAAYZ + DIABZY + DIACZY)
            DADZZ = - FAD*SAZ*SDZ*(DIAAZZ + DIABZZ + DIACZZ)
C
            DBDXX = - FBD*SBX*SDX*(DIABXX + DIBBXX + DIBCXX)
            DBDXY = -     SBX*SDY*(DIABYX + DIBBXY + DIBCXY)
            DBDXZ = -     SBX*SDZ*(DIABZX + DIBBXZ + DIBCXZ)
            DBDYX = -     SBY*SDX*(DIABXY + DIBBXY + DIBCYX)
            DBDYY = - FBD*SBY*SDY*(DIABYY + DIBBYY + DIBCYY)
            DBDYZ = -     SBY*SDZ*(DIABZY + DIBBYZ + DIBCYZ)
            DBDZX = -     SBZ*SDX*(DIABXZ + DIBBXZ + DIBCZX)
            DBDZY = -     SBZ*SDY*(DIABYZ + DIBBYZ + DIBCZY)
            DBDZZ = - FBD*SBZ*SDZ*(DIABZZ + DIBBZZ + DIBCZZ)
C
            DCDXX = - FCD*SCX*SDX*(DIACXX + DIBCXX + DICCXX)
            DCDXY = -     SCX*SDY*(DIACYX + DIBCYX + DICCXY)
            DCDXZ = -     SCX*SDZ*(DIACZX + DIBCZX + DICCXZ)
            DCDYX = -     SCY*SDX*(DIACXY + DIBCXY + DICCXY)
            DCDYY = - FCD*SCY*SDY*(DIACYY + DIBCYY + DICCYY)
            DCDYZ = -     SCY*SDZ*(DIACZY + DIBCZY + DICCYZ)
            DCDZX = -     SCZ*SDX*(DIACXZ + DIBCXZ + DICCXZ)
            DCDZY = -     SCZ*SDY*(DIACYZ + DIBCYZ + DICCYZ)
            DCDZZ = - FCD*SCZ*SDZ*(DIACZZ + DIBCZZ + DICCZZ)
C
            DDDXX =         (DIAAXX + DIABXX + DIACXX + DIABXX + DIBBXX
     *                              + DIBCXX + DIACXX + DIBCXX + DICCXX)
            DDDXY = SDX*SDY*(DIAAXY + DIABXY + DIACXY + DIABYX + DIBBXY
     *                              + DIBCXY + DIACYX + DIBCYX + DICCXY)
            DDDXZ = SDX*SDZ*(DIAAXZ + DIABXZ + DIACXZ + DIABZX + DIBBXZ
     *                              + DIBCXZ + DIACZX + DIBCZX + DICCXZ)
            DDDYY =         (DIAAYY + DIABYY + DIACYY + DIABYY + DIBBYY
     *                              + DIBCYY + DIACYY + DIBCYY + DICCYY)
            DDDYZ = SDY*SDZ*(DIAAYZ + DIABYZ + DIACYZ + DIABZY + DIBBYZ
     *                              + DIBCYZ + DIACZY + DIBCZY + DICCYZ)
            DDDZZ =         (DIAAZZ + DIABZZ + DIACZZ + DIABZZ + DIBBZZ
     *                              + DIBCZZ + DIACZZ + DIBCZZ + DICCZZ)
C
C           Loop over irreps
C
            DO 300 IREP = 0, MAXLOP
               CHIA = PT(IBTAND(ISOA,IREP))
               CHIB = PT(IBTAND(ISOB,IREP))
               CHIC = PT(IBTAND(ISOC,IREP))
               CHID = PT(IBTAND(ISOD,IREP))
               CHIAB = CHIA*CHIB
               CHIAC = CHIA*CHIC
               CHIAD = CHIA*CHID
               CHIBC = CHIB*CHIC
               CHIBD = CHIB*CHID
               CHICD = CHIC*CHID
               IAX  = IPTCNT(NAX,IREP)
               IAY  = IPTCNT(NAY,IREP)
               IAZ  = IPTCNT(NAZ,IREP)
               IBX  = IPTCNT(NBX,IREP)
               IBY  = IPTCNT(NBY,IREP)
               IBZ  = IPTCNT(NBZ,IREP)
               ICX  = IPTCNT(NCX,IREP)
               ICY  = IPTCNT(NCY,IREP)
               ICZ  = IPTCNT(NCZ,IREP)
               IDX  = IPTCNT(NDX,IREP)
               IDY  = IPTCNT(NDY,IREP)
               IDZ  = IPTCNT(NDZ,IREP)
               DOAX = IAX .NE. 0
               DOAY = IAY .NE. 0
               DOAZ = IAZ .NE. 0
               DOBX = IBX .NE. 0
               DOBY = IBY .NE. 0
               DOBZ = IBZ .NE. 0
               DOCX = ICX .NE. 0
               DOCY = ICY .NE. 0
               DOCZ = ICZ .NE. 0
               DODX = IDX .NE. 0
               DODY = IDY .NE. 0
               DODZ = IDZ .NE. 0
C
C              A-A electron-repulsion Hessian elements:
C
               IF (DOAX)
     *            HESSEE(IAX,IAX) = HESSEE(IAX,IAX) + DAAXX
               IF (DOAX.AND.DOAY)
     *            HESSEE(IAX,IAY) = HESSEE(IAX,IAY) + DAAXY
               IF (DOAX.AND.DOAZ)
     *            HESSEE(IAX,IAZ) = HESSEE(IAX,IAZ) + DAAXZ
               IF (DOAY)
     *            HESSEE(IAY,IAY) = HESSEE(IAY,IAY) + DAAYY
               IF (DOAY.AND.DOAZ)
     *            HESSEE(IAY,IAZ) = HESSEE(IAY,IAZ) + DAAYZ
               IF (DOAZ)
     *            HESSEE(IAZ,IAZ) = HESSEE(IAZ,IAZ) + DAAZZ
C
C              A-B electron-repulsion Hessian elements:
C
               IF (DOAX.AND.DOBX)
     *            HESSEE(IAX,IBX) = HESSEE(IAX,IBX) + CHIAB*DABXX
               IF (DOAX.AND.DOBY)
     *            HESSEE(IAX,IBY) = HESSEE(IAX,IBY) + CHIAB*DABXY
               IF (DOAX.AND.DOBZ)
     *            HESSEE(IAX,IBZ) = HESSEE(IAX,IBZ) + CHIAB*DABXZ
               IF (DOAY.AND.DOBX)
     *            HESSEE(IAY,IBX) = HESSEE(IAY,IBX) + CHIAB*DABYX
               IF (DOAY.AND.DOBY)
     *            HESSEE(IAY,IBY) = HESSEE(IAY,IBY) + CHIAB*DABYY
               IF (DOAY.AND.DOBZ)
     *            HESSEE(IAY,IBZ) = HESSEE(IAY,IBZ) + CHIAB*DABYZ
               IF (DOAZ.AND.DOBX)
     *            HESSEE(IAZ,IBX) = HESSEE(IAZ,IBX) + CHIAB*DABZX
               IF (DOAZ.AND.DOBY)
     *            HESSEE(IAZ,IBY) = HESSEE(IAZ,IBY) + CHIAB*DABZY
               IF (DOAZ.AND.DOBZ)
     *            HESSEE(IAZ,IBZ) = HESSEE(IAZ,IBZ) + CHIAB*DABZZ
C
C              A-C electron-repulsion Hessian elements:
C
               IF (DOAX.AND.DOCX)
     *            HESSEE(IAX,ICX) = HESSEE(IAX,ICX) + CHIAC*DACXX
               IF (DOAX.AND.DOCY)
     *            HESSEE(IAX,ICY) = HESSEE(IAX,ICY) + CHIAC*DACXY
               IF (DOAX.AND.DOCZ)
     *            HESSEE(IAX,ICZ) = HESSEE(IAX,ICZ) + CHIAC*DACXZ
               IF (DOAY.AND.DOCX)
     *            HESSEE(IAY,ICX) = HESSEE(IAY,ICX) + CHIAC*DACYX
               IF (DOAY.AND.DOCY)
     *            HESSEE(IAY,ICY) = HESSEE(IAY,ICY) + CHIAC*DACYY
               IF (DOAY.AND.DOCZ)
     *            HESSEE(IAY,ICZ) = HESSEE(IAY,ICZ) + CHIAC*DACYZ
               IF (DOAZ.AND.DOCX)
     *            HESSEE(IAZ,ICX) = HESSEE(IAZ,ICX) + CHIAC*DACZX
               IF (DOAZ.AND.DOCY)
     *            HESSEE(IAZ,ICY) = HESSEE(IAZ,ICY) + CHIAC*DACZY
               IF (DOAZ.AND.DOCZ)
     *            HESSEE(IAZ,ICZ) = HESSEE(IAZ,ICZ) + CHIAC*DACZZ
C
C              A-D electron-repulsion Hessian elements:
C
               IF (DOAX.AND.DODX)
     *            HESSEE(IAX,IDX) = HESSEE(IAX,IDX) + CHIAD*DADXX
               IF (DOAX.AND.DODY)
     *            HESSEE(IAX,IDY) = HESSEE(IAX,IDY) + CHIAD*DADXY
               IF (DOAX.AND.DODZ)
     *            HESSEE(IAX,IDZ) = HESSEE(IAX,IDZ) + CHIAD*DADXZ
               IF (DOAY.AND.DODX)
     *            HESSEE(IAY,IDX) = HESSEE(IAY,IDX) + CHIAD*DADYX
               IF (DOAY.AND.DODY)
     *            HESSEE(IAY,IDY) = HESSEE(IAY,IDY) + CHIAD*DADYY
               IF (DOAY.AND.DODZ)
     *            HESSEE(IAY,IDZ) = HESSEE(IAY,IDZ) + CHIAD*DADYZ
               IF (DOAZ.AND.DODX)
     *            HESSEE(IAZ,IDX) = HESSEE(IAZ,IDX) + CHIAD*DADZX
               IF (DOAZ.AND.DODY)
     *            HESSEE(IAZ,IDY) = HESSEE(IAZ,IDY) + CHIAD*DADZY
               IF (DOAZ.AND.DODZ)
     *            HESSEE(IAZ,IDZ) = HESSEE(IAZ,IDZ) + CHIAD*DADZZ
C
C              B-B electron-repulsion Hessian elements:
C
               IF (DOBX)
     *            HESSEE(IBX,IBX) = HESSEE(IBX,IBX) + DBBXX
               IF (DOBX.AND.DOBY)
     *            HESSEE(IBX,IBY) = HESSEE(IBX,IBY) + DBBXY
               IF (DOBX.AND.DOBZ)
     *            HESSEE(IBX,IBZ) = HESSEE(IBX,IBZ) + DBBXZ
               IF (DOBY)
     *            HESSEE(IBY,IBY) = HESSEE(IBY,IBY) + DBBYY
               IF (DOBY.AND.DOBZ)
     *            HESSEE(IBY,IBZ) = HESSEE(IBY,IBZ) + DBBYZ
               IF (DOBZ)
     *            HESSEE(IBZ,IBZ) = HESSEE(IBZ,IBZ) + DBBZZ
C
C              B-C electron-repulsion Hessian elements:
C
               IF (DOBX.AND.DOCX)
     *            HESSEE(IBX,ICX) = HESSEE(IBX,ICX) + CHIBC*DBCXX
               IF (DOBX.AND.DOCY)
     *            HESSEE(IBX,ICY) = HESSEE(IBX,ICY) + CHIBC*DBCXY
               IF (DOBX.AND.DOCZ)
     *            HESSEE(IBX,ICZ) = HESSEE(IBX,ICZ) + CHIBC*DBCXZ
               IF (DOBY.AND.DOCX)
     *            HESSEE(IBY,ICX) = HESSEE(IBY,ICX) + CHIBC*DBCYX
               IF (DOBY.AND.DOCY)
     *            HESSEE(IBY,ICY) = HESSEE(IBY,ICY) + CHIBC*DBCYY
               IF (DOBY.AND.DOCZ)
     *            HESSEE(IBY,ICZ) = HESSEE(IBY,ICZ) + CHIBC*DBCYZ
               IF (DOBZ.AND.DOCX)
     *            HESSEE(IBZ,ICX) = HESSEE(IBZ,ICX) + CHIBC*DBCZX
               IF (DOBZ.AND.DOCY)
     *            HESSEE(IBZ,ICY) = HESSEE(IBZ,ICY) + CHIBC*DBCZY
               IF (DOBZ.AND.DOCZ)
     *            HESSEE(IBZ,ICZ) = HESSEE(IBZ,ICZ) + CHIBC*DBCZZ
C
C              B-D electron-repulsion Hessian elements:
C
               IF (DOBX.AND.DODX)
     *            HESSEE(IBX,IDX) = HESSEE(IBX,IDX) + CHIBD*DBDXX
               IF (DOBX.AND.DODY)
     *            HESSEE(IBX,IDY) = HESSEE(IBX,IDY) + CHIBD*DBDXY
               IF (DOBX.AND.DODZ)
     *            HESSEE(IBX,IDZ) = HESSEE(IBX,IDZ) + CHIBD*DBDXZ
               IF (DOBY.AND.DODX)
     *            HESSEE(IBY,IDX) = HESSEE(IBY,IDX) + CHIBD*DBDYX
               IF (DOBY.AND.DODY)
     *            HESSEE(IBY,IDY) = HESSEE(IBY,IDY) + CHIBD*DBDYY
               IF (DOBY.AND.DODZ)
     *            HESSEE(IBY,IDZ) = HESSEE(IBY,IDZ) + CHIBD*DBDYZ
               IF (DOBZ.AND.DODX)
     *            HESSEE(IBZ,IDX) = HESSEE(IBZ,IDX) + CHIBD*DBDZX
               IF (DOBZ.AND.DODY)
     *            HESSEE(IBZ,IDY) = HESSEE(IBZ,IDY) + CHIBD*DBDZY
               IF (DOBZ.AND.DODZ)
     *            HESSEE(IBZ,IDZ) = HESSEE(IBZ,IDZ) + CHIBD*DBDZZ
C
C              C-C electron-repulsion Hessian elements:
C
               IF (DOCX)
     *            HESSEE(ICX,ICX) = HESSEE(ICX,ICX) + DCCXX
               IF (DOCX.AND.DOCY)
     *            HESSEE(ICX,ICY) = HESSEE(ICX,ICY) + DCCXY
               IF (DOCX.AND.DOCZ)
     *            HESSEE(ICX,ICZ) = HESSEE(ICX,ICZ) + DCCXZ
               IF (DOCY)
     *            HESSEE(ICY,ICY) = HESSEE(ICY,ICY) + DCCYY
               IF (DOCY.AND.DOCZ)
     *            HESSEE(ICY,ICZ) = HESSEE(ICY,ICZ) + DCCYZ
               IF (DOCZ)
     *            HESSEE(ICZ,ICZ) = HESSEE(ICZ,ICZ) + DCCZZ
C
C              C-D electron-repulsion Hessian elements:
C
               IF (DOCX.AND.DODX)
     *            HESSEE(ICX,IDX) = HESSEE(ICX,IDX) + CHICD*DCDXX
               IF (DOCX.AND.DODY)
     *            HESSEE(ICX,IDY) = HESSEE(ICX,IDY) + CHICD*DCDXY
               IF (DOCX.AND.DODZ)
     *            HESSEE(ICX,IDZ) = HESSEE(ICX,IDZ) + CHICD*DCDXZ
               IF (DOCY.AND.DODX)
     *            HESSEE(ICY,IDX) = HESSEE(ICY,IDX) + CHICD*DCDYX
               IF (DOCY.AND.DODY)
     *            HESSEE(ICY,IDY) = HESSEE(ICY,IDY) + CHICD*DCDYY
               IF (DOCY.AND.DODZ)
     *            HESSEE(ICY,IDZ) = HESSEE(ICY,IDZ) + CHICD*DCDYZ
               IF (DOCZ.AND.DODX)
     *            HESSEE(ICZ,IDX) = HESSEE(ICZ,IDX) + CHICD*DCDZX
               IF (DOCZ.AND.DODY)
     *            HESSEE(ICZ,IDY) = HESSEE(ICZ,IDY) + CHICD*DCDZY
               IF (DOCZ.AND.DODZ)
     *            HESSEE(ICZ,IDZ) = HESSEE(ICZ,IDZ) + CHICD*DCDZZ
C
C              D-D electron-repulsion Hessian elements:
C
               IF (DODX)
     *            HESSEE(IDX,IDX) = HESSEE(IDX,IDX) + DDDXX
               IF (DODX.AND.DODY)
     *            HESSEE(IDX,IDY) = HESSEE(IDX,IDY) + DDDXY
               IF (DODX.AND.DODZ)
     *            HESSEE(IDX,IDZ) = HESSEE(IDX,IDZ) + DDDXZ
               IF (DODY)
     *            HESSEE(IDY,IDY) = HESSEE(IDY,IDY) + DDDYY
               IF (DODY.AND.DODZ)
     *            HESSEE(IDY,IDZ) = HESSEE(IDY,IDZ) + DDDYZ
               IF (DODZ)
     *            HESSEE(IDZ,IDZ) = HESSEE(IDZ,IDZ) + DDDZZ
  300       CONTINUE
         END IF
      END IF
      IF (IPRINT .GT. 05) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 1110) TWOCEN, THRCEN
         WRITE (LUPRI, 1120) NCENT1, NCENT2, NCENT3, NCENT4
         WRITE (LUPRI,'(A,4I5)') ' ISOA/B/C/D ',ISOA,ISOB,ISOC,ISOD
         WRITE (LUPRI,'(A,3I5)') ' IAX/Y/Z ',IAX,IAY,IAZ
         WRITE (LUPRI,'(A,3I5)') ' IBX/Y/Z ',IBX,IBY,IBZ
         WRITE (LUPRI,'(A,3I5)') ' ICX/Y/Z ',ICX,ICY,ICZ
         WRITE (LUPRI,'(A,3I5)') ' IDX/Y/Z ',IDX,IDY,IDZ
         WRITE (LUPRI,'(A,3F12.6)') ' SAX/Y/Z ',SAX,SAY,SAZ
         WRITE (LUPRI,'(A,3F12.6)') ' SBX/Y/Z ',SBX,SBY,SBZ
         WRITE (LUPRI,'(A,3F12.6)') ' SCX/Y/Z ',SCX,SCY,SCZ
         WRITE (LUPRI,'(A,3F12.6)') ' SDX/Y/Z ',SDX,SDY,SDZ
         WRITE (LUPRI, 4000) (DERIV(I),I = 1, NINTYP)
         NCDEP3 = 3*NUCDEP
         IF (IPRINT .GE. 10) THEN
            IF (DERONE) THEN
               WRITE(LUPRI, 4010)
               WRITE(LUPRI, 4020) (GRADEE(I),I=1,NCDEP3)
            END IF
            IF (DERTWO) THEN
               WRITE(LUPRI, 4030)
               DO 800 I = 1,NCDEP3
                  WRITE (LUPRI, 4020)
     *                     (HESSEE(I,J) + HESSEE(J,I), J = 1,I)
  800          CONTINUE
            END IF
         END IF
      END IF
      RETURN
 1000 FORMAT (//,1X,' <<<<<<<<<< Subroutine DEROUT >>>>>>>>',/)
 1110 FORMAT (2X,'TWOCEN,...:  ',2L5)
 1120 FORMAT (2X,'NCENT1/2/3/4:',4I7)
 4000 FORMAT (2X,'DERIV        ',4E15.6,/(15X,4E15.6))
 4010 FORMAT (//,1X,' Two-electron integral gradient ',/)
 4020 FORMAT (6F12.6)
 4030 FORMAT (//,1X,' Two-electron integral Hessian ',/)
      END
