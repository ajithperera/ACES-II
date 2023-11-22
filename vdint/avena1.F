      SUBROUTINE AVENA1(KHKTA,KHKTB,IDENA,IDENB,ICENTA,NCENTA,
     &                  NATOMC,SECOND,LDIAG,WORK1,ENERNA,
     &                  GRADNA,HESSNA,NCOORD,ECPOPX,ECPOPY,ECPOPZ)
C
C  CALCULATES THE VARIOUS EXPECTATION VALUES FOR THE NUCLEAR
C  ATTRACTION INTEGRALS FOR ONE-CENTER INTEGRALS
C
CEND
C
C  ADAPTED TO THE CRAPS ENVIRONMENT  OCT/90/JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
#include "mxcent.par"
#include "baslims.par"
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN*
     *           MXCONT*MXCONT)
      PARAMETER (D0 = 0.00 D00, D1 = 1.00 D00, D2 = 2.00D00)
      LOGICAL LDIAG
      DIMENSION WORK1(1),GRADNA(NCOORD),HESSNA(NCOORD,NCOORD)
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
      COMMON /CENTC/ SIGNCX(MXCENT), SIGNCY(MXCENT), SIGNCZ(MXCENT),
     *               NCENTC(MXCENT), JSYMC(MXCENT),  JCENTC(MXCENT),
     *               ICXVEC(MXCENT), ICYVEC(MXCENT), ICZVEC(MXCENT)
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
      COMMON/GENCON/NRCA,NRCB,CONTA(MXCONT),CONTB(MXCONT)
      LOGICAL SECOND, DOAX, DOAY, DOAZ, DOCX, DOCY, DOCZ
      COMMON/DF/ DSHELL(MXAQNS),FSHELL(MXAQNS)
      COMMON /ECP   / IECP
CCH---------------------------------------------------------------------
C     ECPCONTRIBUTION TO GRAD
CCH---------------------------------------------------------------------
      DIMENSION ECPOPX(KHKTB*KHKTA*NATOMC*MXCONT*MXCONT)
      DIMENSION ECPOPY(KHKTB*KHKTA*NATOMC*MXCONT*MXCONT)
      DIMENSION ECPOPZ(KHKTB*KHKTA*NATOMC*MXCONT*MXCONT)
CCH---------------------------------------------------------------------
C
      DATA HALF,TWO /0.5D0,2.D0/
C
      IBTAND(I,J) = AND(I,J)
      ITRI(I,J) = MAX(I,J)*(MAX(I,J) - 1)/2 + MIN(I,J)
C
C     Get density elements
C
      IOFF=0
      ISKIP=NRCA*NRCB
      IF(LDIAG) ISKIP=NRCA*(NRCA+1)/2
      DO 200 IRCA=1,NRCA
       MAXB=NRCB
       IF(LDIAG) MAXB=IRCA
       DO 200 IRCB=1,MAXB
        IOFF=IOFF+1
        ICOMP = IOFF
        DO 100 IORBA = IDENA + IRCA, IDENA + NRCA*KHKTA,NRCA
        DO 150 IORBB = IDENB + IRCB, IDENB + NRCB*KHKTB,NRCB
         DSHELL(ICOMP) = WORK1(KDEN-1+ITRI(IORBA,IORBB))
         ICOMP = ICOMP + ISKIP
150     CONTINUE
100    CONTINUE
       IF(LDIAG.AND.(IRCA.EQ.IRCB)) THEN
cch     multiply ALL elements by 1/2
        CALL SSCAL(KHKTA*KHKTB,HALF,DSHELL(IOFF),ISKIP)
cch     correction: multiply diagonal elements by 2
        CALL SSCAL(KHKTA,TWO,DSHELL(IOFF),ISKIP*(KHKTA+1))
       ENDIF
200   CONTINUE
      MAXCMP = ISKIP*KHKTA*KHKTB
C
      NAX = 3*NCENTA - 2
      NAY = 3*NCENTA - 1
      NAZ = 3*NCENTA
C
C     Run over atoms in the operator
C
#ifdef _OLDECP
CCH---------------------------------------------------------------------
C     ADD ECPCONTRIBUTIONS
CCH---------------------------------------------------------------------
      IF (IECP.EQ.1) THEN
        IOFF=KHKTB*KHKTA*NATOMC
        IWODIM=NATOMC*NRCA*NRCB
        IF (LDIAG) IWODIM=NATOMC*NRCA*(NRCA+1)/2
        ICOUNT=0
        DO 5000 IRCA=1,NRCA
          MAXB=NRCB
          IF(LDIAG) MAXB=IRCA
          DO 5100 IRCB=1,MAXB
            IPOINT=IOFF*ICOUNT
            ICOUNT=ICOUNT+1
            CALL ECPADD(WORK1(IA000X+1),ECPOPX(IPOINT+1),KHKTA,KHKTB,
     &                  NATOMC,NRCB,NRCA,IWODIM,ICOUNT,LDIAG)
            CALL ECPADD(WORK1(IA000Y+1),ECPOPY(IPOINT+1),KHKTA,KHKTB,
     &                  NATOMC,NRCB,NRCA,IWODIM,ICOUNT,LDIAG)
            CALL ECPADD(WORK1(IA000Z+1),ECPOPZ(IPOINT+1),KHKTA,KHKTB,
     &                  NATOMC,NRCB,NRCA,IWODIM,ICOUNT,LDIAG)
 5100     CONTINUE
 5000   CONTINUE
      ENDIF
CCH---------------------------------------------------------------------
#endif
      DO 250 IATOM = 1,NATOMC
       ICENTC = NCENTC(IATOM)
       IF (ICENTC .NE. ICENTA) THEN
        DA0000 = SDOT(MAXCMP,DSHELL,1,WORK1(IA0000+IATOM),NATOMC)
        DA000X = SDOT(MAXCMP,DSHELL,1,WORK1(IA000X+IATOM),NATOMC)
        DA000Y = SDOT(MAXCMP,DSHELL,1,WORK1(IA000Y+IATOM),NATOMC)
        DA000Z = SDOT(MAXCMP,DSHELL,1,WORK1(IA000Z+IATOM),NATOMC)
C
        IF(SECOND) THEN
         DA00XX = SDOT(MAXCMP,DSHELL,1,WORK1(IA00XX+IATOM),NATOMC)
         DA00XY = SDOT(MAXCMP,DSHELL,1,WORK1(IA00XY+IATOM),NATOMC)
         DA00XZ = SDOT(MAXCMP,DSHELL,1,WORK1(IA00XZ+IATOM),NATOMC)
         DA00YY = SDOT(MAXCMP,DSHELL,1,WORK1(IA00YY+IATOM),NATOMC)
         DA00YZ = SDOT(MAXCMP,DSHELL,1,WORK1(IA00YZ+IATOM),NATOMC)
         DA00ZZ = SDOT(MAXCMP,DSHELL,1,WORK1(IA00ZZ+IATOM),NATOMC)
        ENDIF
C
        KCENTC = JCENTC(IATOM)
        NCX    = 3*KCENTC - 2
        NCY    = 3*KCENTC - 1
        NCZ    = 3*KCENTC
        SCX    = SIGNCX(IATOM)
        SCY    = SIGNCY(IATOM)
        SCZ    = SIGNCZ(IATOM)
C
C  Energy 
C
        ENERNA = ENERNA + DA0000
C
C  Gradient 
C
       IAX  = IPTCNT(NAX,0)
       IAY  = IPTCNT(NAY,0)
       IAZ  = IPTCNT(NAZ,0)
       ICX  = IPTCNT(NCX,0)
       ICY  = IPTCNT(NCY,0)
       ICZ  = IPTCNT(NCZ,0)
       DOAX = IAX .NE. 0
       DOAY = IAY .NE. 0
       DOAZ = IAZ .NE. 0
       DOCX = ICX .NE. 0
       DOCY = ICY .NE. 0
       DOCZ = ICZ .NE. 0
C
C  A nuclear-attraction gradient elements:
C
       IF (DOAX) GRADNA(IAX) = GRADNA(IAX) - DA000X
       IF (DOAY) GRADNA(IAY) = GRADNA(IAY) - DA000Y
       IF (DOAZ) GRADNA(IAZ) = GRADNA(IAZ) - DA000Z
C
C   C nuclear-attraction gradient elements:
C
       IF (DOCX) GRADNA(ICX) = GRADNA(ICX) + SCX*DA000X
       IF (DOCY) GRADNA(ICY) = GRADNA(ICY) + SCY*DA000Y
       IF (DOCZ) GRADNA(ICZ) = GRADNA(ICZ) + SCZ*DA000Z
C
C  Hessian *****
C
            IF (SECOND) THEN
               ISYMPC = JSYMC(IATOM)
               FAC = D1
               IF (NCENTA .EQ. KCENTC) FAC = D2
C
C              Run over irreps
C
               DO 400 IREP = 0, MAXLOP
                  CHIC = PT(IBTAND(ISYMPC,IREP))
                  CSCX = CHIC*SCX
                  CSCY = CHIC*SCY
                  CSCZ = CHIC*SCZ
                  IAX  = IPTCNT(NAX,IREP)
                  IAY  = IPTCNT(NAY,IREP)
                  IAZ  = IPTCNT(NAZ,IREP)
                  ICX  = IPTCNT(NCX,IREP)
                  ICY  = IPTCNT(NCY,IREP)
                  ICZ  = IPTCNT(NCZ,IREP)
                  DOAX = IAX .NE. 0
                  DOAY = IAY .NE. 0
                  DOAZ = IAZ .NE. 0
                  DOCX = ICX .NE. 0
                  DOCY = ICY .NE. 0
                  DOCZ = ICZ .NE. 0
C
C                 A-A nuclear-attraction Hessian elements:
C
                  IF (DOAX)          HESSNA(IAX,IAX) =
     *                               HESSNA(IAX,IAX) + DA00XX
                  IF (DOAX.AND.DOAY) HESSNA(IAX,IAY) =
     *                               HESSNA(IAX,IAY) + DA00XY
                  IF (DOAX.AND.DOAZ) HESSNA(IAX,IAZ) =
     *                               HESSNA(IAX,IAZ) + DA00XZ
                  IF (DOAY)          HESSNA(IAY,IAY) =
     *                               HESSNA(IAY,IAY) + DA00YY
                  IF (DOAY.AND.DOAZ) HESSNA(IAY,IAZ) =
     *                               HESSNA(IAY,IAZ) + DA00YZ
                  IF (DOAZ)          HESSNA(IAZ,IAZ) =
     *                               HESSNA(IAZ,IAZ) + DA00ZZ
C
C                 A-C nuclear-attraction Hessian elements:
C
                  IF (DOAX.AND.DOCX) HESSNA(IAX,ICX) =
     *                               HESSNA(IAX,ICX) - FAC*CSCX*DA00XX
                  IF (DOAX.AND.DOCY) HESSNA(IAX,ICY) =
     *                               HESSNA(IAX,ICY) - CSCY*DA00XY
                  IF (DOAX.AND.DOCZ) HESSNA(IAX,ICZ) =
     *                               HESSNA(IAX,ICZ) - CSCZ*DA00XZ
                  IF (DOAY.AND.DOCX) HESSNA(IAY,ICX) =
     *                               HESSNA(IAY,ICX) - CSCX*DA00XY
                  IF (DOAY.AND.DOCY) HESSNA(IAY,ICY) =
     *                               HESSNA(IAY,ICY) - FAC*CSCY*DA00YY
                  IF (DOAY.AND.DOCZ) HESSNA(IAY,ICZ) =
     *                               HESSNA(IAY,ICZ) - CSCZ*DA00YZ
                  IF (DOAZ.AND.DOCX) HESSNA(IAZ,ICX) =
     *                               HESSNA(IAZ,ICX) - CSCX*DA00XZ
                  IF (DOAZ.AND.DOCY) HESSNA(IAZ,ICY) =
     *                               HESSNA(IAZ,ICY) - CSCY*DA00YZ
                  IF (DOAZ.AND.DOCZ) HESSNA(IAZ,ICZ) =
     *                               HESSNA(IAZ,ICZ) - FAC*CSCZ*DA00ZZ
C
C                 C-C nuclear-attraction Hessian elements:
C
                  IF (DOCX)          HESSNA(ICX,ICX) =
     *                               HESSNA(ICX,ICX) + DA00XX
                  IF (DOCX.AND.DOCY) HESSNA(ICX,ICY) =
     *                               HESSNA(ICX,ICY) + SCX*SCY*DA00XY
                  IF (DOCX.AND.DOCZ) HESSNA(ICX,ICZ) =
     *                               HESSNA(ICX,ICZ) + SCX*SCZ*DA00XZ
                  IF (DOCY)          HESSNA(ICY,ICY) =
     *                               HESSNA(ICY,ICY) + DA00YY
                  IF (DOCY.AND.DOCZ) HESSNA(ICY,ICZ) =
     *                               HESSNA(ICY,ICZ) + SCY*SCZ*DA00YZ
                  IF (DOCZ)          HESSNA(ICZ,ICZ) =
     *                               HESSNA(ICZ,ICZ) + DA00ZZ
  400          CONTINUE
            END IF
         END IF
  250 CONTINUE
      RETURN
      END
