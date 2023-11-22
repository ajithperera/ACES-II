      SUBROUTINE AVENA2(KHKTA,KHKTB,IDENA,IDENB,ICENTA,ICENTB,
     &                  NCENTA,NCENTB,SBX,SBY,SBZ,
     &                  NATOMC,ISYMOP,SECOND,LDIAG,WORK1,
     &                  ENERNA,GRADNA,HESSNA,NCOORD,
     &                  ECPOPX,ECPOPY,ECPOPZ,ECPBFX,ECPBFY,ECPBFZ)
C
CEND
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
#include "mxcent.par"
#include "baslims.par"
      PARAMETER (D0 = 0.00 D00, D1 = 1.00 D00, D2 = 2.00D00)
      LOGICAL LDIAG
      DIMENSION WORK1(1)
      DIMENSION GRADNA(NCOORD),HESSNA(NCOORD,NCOORD)

      COMMON /CENTC/ SIGNCX(MXCENT), SIGNCY(MXCENT), SIGNCZ(MXCENT),
     *               NCENTC(MXCENT), JSYMC(MXCENT),  JCENTC(MXCENT),
     *               ICXVEC(MXCENT), ICYVEC(MXCENT), ICZVEC(MXCENT)
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN
     *           *MXCONT*MXCONT)
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
      COMMON /ADER/ ADER0 (MXAQNS)
      COMMON/POINTER/ IS0000, IS000X, IS000Y, IS000Z,
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
      LOGICAL SECOND, DOAX, DOAY, DOAZ, DOBX, DOBY, DOBZ,
     *        DOCX, DOCY, DOCZ
      COMMON/DF/DSHELL(MXAQNS),FSHELL(MXAQNS)
      COMMON /ECP   / IECP
CCH---------------------------------------------------------------------
C     ECPCONTRIBUTION TO GRAD
CCH---------------------------------------------------------------------
      DIMENSION ECPOPX(KHKTB*KHKTA*NATOMC*MXCONT*MXCONT)
      DIMENSION ECPOPY(KHKTB*KHKTA*NATOMC*MXCONT*MXCONT)
      DIMENSION ECPOPZ(KHKTB*KHKTA*NATOMC*MXCONT*MXCONT)
      DIMENSION ECPBFX(KHKTB*KHKTA*NATOMC*MXCONT*MXCONT)
      DIMENSION ECPBFY(KHKTB*KHKTA*NATOMC*MXCONT*MXCONT)
      DIMENSION ECPBFZ(KHKTB*KHKTA*NATOMC*MXCONT*MXCONT)
CCH---------------------------------------------------------------------

C
      DATA HALF,TWO /0.5D0,2.D0/
c
      IBTAND(I,J) = AND(I,J)
      ITRI(I,J) = MAX(I,J)*(MAX(I,J) - 1)/2 + MIN(I,J)
C
C     Get density elements
C
      IOFF=0
      ISKIP=NRCA*NRCB
      IF(LDIAG) ISKIP=NRCA*(NRCA+1)/2
      MAXCMP = ISKIP*KHKTA*KHKTB
C
      NAX = 3*NCENTA - 2
      NAY = 3*NCENTA - 1
      NAZ = 3*NCENTA
      NBX = 3*NCENTB - 2
      NBY = 3*NCENTB - 1
      NBZ = 3*NCENTB
      FAB = D1
      IF (NCENTA .EQ. NCENTB) FAB = D2
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
     &                        NATOMC,NRCB,NRCA,IWODIM,ICOUNT,LDIAG)
            CALL ECPADD(WORK1(IA000Y+1),ECPOPY(IPOINT+1),KHKTA,KHKTB,
     &                        NATOMC,NRCB,NRCA,IWODIM,ICOUNT,LDIAG)
            CALL ECPADD(WORK1(IA000Z+1),ECPOPZ(IPOINT+1),KHKTA,KHKTB,
     &                        NATOMC,NRCB,NRCA,IWODIM,ICOUNT,LDIAG)
            CALL ECPADD(WORK1(IA0X00+1),ECPBFX(IPOINT+1),KHKTA,KHKTB,
     &                        NATOMC,NRCB,NRCA,IWODIM,ICOUNT,LDIAG)
            CALL ECPADD(WORK1(IA0Y00+1),ECPBFY(IPOINT+1),KHKTA,KHKTB,
     &                        NATOMC,NRCB,NRCA,IWODIM,ICOUNT,LDIAG)
            CALL ECPADD(WORK1(IA0Z00+1),ECPBFZ(IPOINT+1),KHKTA,KHKTB,
     &                        NATOMC,NRCB,NRCA,IWODIM,ICOUNT,LDIAG)
 5100     CONTINUE
 5000   CONTINUE
      ENDIF
cchdebug      
c      write(*,*) 'avena2'
c      write(*,*) 'basisfunct.'
c      DO 5001 IRCA=1,NRCA
c        MAXB=NRCB
c        IF(LDIAG) MAXB=IRCA
c        DO 5101 IRCB=1,MAXB
c          IPOINT=IOFF*ICOUNT
c          ICOUNT=ICOUNT+1
c          CALL ECPADD(WORK1(IA0X00+1),ECPBFX(IPOINT+1),KHKTA,KHKTB,
c     &                      NATOMC,NRCB,NRCA,IWODIM,ICOUNT,LDIAG)
c          CALL ECPADD(WORK1(IA0Y00+1),ECPBFY(IPOINT+1),KHKTA,KHKTB,
c     &                      NATOMC,NRCB,NRCA,IWODIM,ICOUNT,LDIAG)
c          CALL ECPADD(WORK1(IA0Z00+1),ECPBFZ(IPOINT+1),KHKTA,KHKTB,
c     &                      NATOMC,NRCB,NRCA,IWODIM,ICOUNT,LDIAG)
c 5101   CONTINUE
c 5001 CONTINUE
cc      write(*,*) 'operator'
c      ICOUNT=0
c      DO 5000 IRCA=1,NRCA
c        MAXB=NRCB
c        IF(LDIAG) MAXB=IRCA
c        DO 5100 IRCB=1,MAXB
c          IPOINT=IOFF*ICOUNT
c          ICOUNT=ICOUNT+1
c          CALL ECPADD(WORK1(IA000X+1),ECPOPX(IPOINT+1),KHKTA,KHKTB,
c     &                      NATOMC,NRCB,NRCA,IWODIM,ICOUNT,LDIAG)
c          CALL ECPADD(WORK1(IA000Y+1),ECPOPY(IPOINT+1),KHKTA,KHKTB,
c     &                      NATOMC,NRCB,NRCA,IWODIM,ICOUNT,LDIAG)
c          CALL ECPADD(WORK1(IA000Z+1),ECPOPZ(IPOINT+1),KHKTA,KHKTB,
c     &                      NATOMC,NRCB,NRCA,IWODIM,ICOUNT,LDIAG)
c 5100   CONTINUE
c 5000 CONTINUE

c      do 99 iatom=1,natomc
c      write(*,*)
c      write(*,*) 'Atomnr.',iatom
c      write(*,*)
c      write(*,*) 'Type ia0x000'
c      write(*,*)
c      write(*,*) 'nucl. in avena2 d/dx',(work1(ia0x00+iatom-1+i),
c     &            i=1,maxcmp*natomc,natomc)
c      write(*,*) 'nucl. in avena2 d/dy',(work1(ia0y00+iatom-1+i),
c     &            i=1,maxcmp*natomc,natomc)
c      write(*,*) 'nucl. in avena2 d/dz',(work1(ia0z00+iatom-1+i),
c     &            i=1,maxcmp*natomc,natomc)
c      write(*,*)
c      write(*,*)
c      write(*,*) 'Type ia0000x'
c      write(*,*)
c      write(*,*) 'nucl. in avena2 d/dx',(work1(ia000x+iatom-1+i),
c     &            i=1,maxcmp*natomc,natomc)
c      write(*,*) 'nucl. in avena2 d/dy',(work1(ia000y+iatom-1+i),
c     &            i=1,maxcmp*natomc,natomc)
c      write(*,*) 'nucl. in avena2 d/dz',(work1(ia000z+iatom-1+i),
c     &            i=1,maxcmp*natomc,natomc)
c      write(*,*) '----------------------------------------------'
c      write(*,*)
c  99  continue
CCH---------------------------------------------------------------------
#endif
C
C    Loop over atoms in operator
C
      DO 250 IATOM = 1,NATOMC
       ICENTC = NCENTC(IATOM)
       ISYMPC = JSYMC(IATOM)
       KCENTC = JCENTC(IATOM)
       NCX    = 3*KCENTC - 2
       NCY    = 3*KCENTC - 1
       NCZ    = 3*KCENTC
       SCX    = SIGNCX(IATOM)
       SCY    = SIGNCY(IATOM)
       SCZ    = SIGNCZ(IATOM)
C
C        Multiply densities and integrals
C
c            IF (SECOND) THEN
c               DAXX00 = DAXX00 + DENS*WORK1(IAXX00 + INT)
c               DAXY00 = DAXY00 + DENS*WORK1(IAXY00 + INT)
c               DAXZ00 = DAXZ00 + DENS*WORK1(IAXZ00 + INT)
c               DAYY00 = DAYY00 + DENS*WORK1(IAYY00 + INT)
c               DAYZ00 = DAYZ00 + DENS*WORK1(IAYZ00 + INT)
c               DAZZ00 = DAZZ00 + DENS*WORK1(IAZZ00 + INT)
c               DA00XX = DA00XX + DENS*WORK1(IA00XX + INT)
c               DA00XY = DA00XY + DENS*WORK1(IA00XY + INT)
c               DA00XZ = DA00XZ + DENS*WORK1(IA00XZ + INT)
c               DA00YY = DA00YY + DENS*WORK1(IA00YY + INT)
c               DA00YZ = DA00YZ + DENS*WORK1(IA00YZ + INT)
c               DA00ZZ = DA00ZZ + DENS*WORK1(IA00ZZ + INT)
c               DA0X0X = DA0X0X + DENS*WORK1(IA0X0X + INT)
c               DA0X0Y = DA0X0Y + DENS*WORK1(IA0X0Y + INT)
c               DA0X0Z = DA0X0Z + DENS*WORK1(IA0X0Z + INT)
c               DA0Y0X = DA0Y0X + DENS*WORK1(IA0Y0X + INT)
c               DA0Y0Y = DA0Y0Y + DENS*WORK1(IA0Y0Y + INT)
c               DA0Y0Z = DA0Y0Z + DENS*WORK1(IA0Y0Z + INT)
c               DA0Z0X = DA0Z0X + DENS*WORK1(IA0Z0X + INT)
c               DA0Z0Y = DA0Z0Y + DENS*WORK1(IA0Z0Y + INT)
c               DA0Z0Z = DA0Z0Z + DENS*WORK1(IA0Z0Z + INT)
c            END IF
c            INT = INT + NATOMC
c  300    CONTINUE
         DA0000=SDOT(MAXCMP,DSHELL,1,WORK1(IA0000+IATOM),NATOMC)
         DA0X00=SDOT(MAXCMP,DSHELL,1,WORK1(IA0X00+IATOM),NATOMC)
         DA0Y00=SDOT(MAXCMP,DSHELL,1,WORK1(IA0Y00+IATOM),NATOMC)
         DA0Z00=SDOT(MAXCMP,DSHELL,1,WORK1(IA0Z00+IATOM),NATOMC)
C
         DA000X=SDOT(MAXCMP,DSHELL,1,WORK1(IA000X+IATOM),NATOMC)
         DA000Y=SDOT(MAXCMP,DSHELL,1,WORK1(IA000Y+IATOM),NATOMC)
         DA000Z=SDOT(MAXCMP,DSHELL,1,WORK1(IA000Z+IATOM),NATOMC)
         IF(SECOND) THEN       
          DAXX00=SDOT(MAXCMP,DSHELL,1,WORK1(IAXX00+IATOM),NATOMC)
          DAXY00=SDOT(MAXCMP,DSHELL,1,WORK1(IAXY00+IATOM),NATOMC)
          DAXZ00=SDOT(MAXCMP,DSHELL,1,WORK1(IAXZ00+IATOM),NATOMC)
          DAYY00=SDOT(MAXCMP,DSHELL,1,WORK1(IAYY00+IATOM),NATOMC)
          DAYZ00=SDOT(MAXCMP,DSHELL,1,WORK1(IAYZ00+IATOM),NATOMC)
          DAZZ00=SDOT(MAXCMP,DSHELL,1,WORK1(IAZZ00+IATOM),NATOMC)
          DA00XX=SDOT(MAXCMP,DSHELL,1,WORK1(IA00XX+IATOM),NATOMC)
          DA00XY=SDOT(MAXCMP,DSHELL,1,WORK1(IA00XY+IATOM),NATOMC)
          DA00XZ=SDOT(MAXCMP,DSHELL,1,WORK1(IA00XZ+IATOM),NATOMC)
          DA00YY=SDOT(MAXCMP,DSHELL,1,WORK1(IA00YY+IATOM),NATOMC)
          DA00YZ=SDOT(MAXCMP,DSHELL,1,WORK1(IA00YZ+IATOM),NATOMC)
          DA00ZZ=SDOT(MAXCMP,DSHELL,1,WORK1(IA00ZZ+IATOM),NATOMC)
          DA0X0X=SDOT(MAXCMP,DSHELL,1,WORK1(IA0X0X+IATOM),NATOMC)
          DA0X0Y=SDOT(MAXCMP,DSHELL,1,WORK1(IA0X0Y+IATOM),NATOMC)
          DA0X0Z=SDOT(MAXCMP,DSHELL,1,WORK1(IA0X0Z+IATOM),NATOMC)
          DA0Y0X=SDOT(MAXCMP,DSHELL,1,WORK1(IA0Y0X+IATOM),NATOMC)
          DA0Y0Y=SDOT(MAXCMP,DSHELL,1,WORK1(IA0Y0Y+IATOM),NATOMC)
          DA0Y0Z=SDOT(MAXCMP,DSHELL,1,WORK1(IA0Y0Z+IATOM),NATOMC)
          DA0Z0X=SDOT(MAXCMP,DSHELL,1,WORK1(IA0Z0X+IATOM),NATOMC)
          DA0Z0Y=SDOT(MAXCMP,DSHELL,1,WORK1(IA0Z0Y+IATOM),NATOMC)
          DA0Z0Z=SDOT(MAXCMP,DSHELL,1,WORK1(IA0Z0Z+IATOM),NATOMC)
         ENDIF
C
C  Undifferentiated Nuclear Attraction Energy
C
         ENERNA = ENERNA + DA0000
C
C  two-center case: C = A *****
C
         IF (ICENTC .EQ. ICENTA) THEN
          TERMX  = DA0X00 + SCX*DA000X
          TERMY  = DA0Y00 + SCY*DA000Y
          TERMZ  = DA0Z00 + SCZ*DA000Z
          IAX    = IPTCNT(NAX,0)
          IAY    = IPTCNT(NAY,0)
          IAZ    = IPTCNT(NAZ,0)
          IBX    = IPTCNT(NBX,0)
          IBY    = IPTCNT(NBY,0)
          IBZ    = IPTCNT(NBZ,0)
          DOAX   = IAX .NE. 0
          DOAY   = IAY .NE. 0
          DOAZ   = IAZ .NE. 0
          DOBX   = IBX .NE. 0
          DOBY   = IBY .NE. 0
          DOBZ   = IBZ .NE. 0
C
C   A nuclear-attraction gradient elements
C
          IF (DOAX) GRADNA(IAX) = GRADNA(IAX) + TERMX
          IF (DOAY) GRADNA(IAY) = GRADNA(IAY) + TERMY
          IF (DOAZ) GRADNA(IAZ) = GRADNA(IAZ) + TERMZ
C
C   B nuclear-attraction gradient elements
C
          IF (DOBX) GRADNA(IBX) = GRADNA(IBX) - SBX*TERMX
          IF (DOBY) GRADNA(IBY) = GRADNA(IBY) - SBY*TERMY
          IF (DOBZ) GRADNA(IBZ) = GRADNA(IBZ) - SBZ*TERMZ
C
C  Second Derivatives:
C
          IF (SECOND) THEN
           TERMXX = DAXX00 + DA0X0X + DA0X0X + DA00XX
           TERMXY = DAXY00 + DA0X0Y + DA0Y0X + DA00XY
           TERMXZ = DAXZ00 + DA0X0Z + DA0Z0X + DA00XZ
           TERMYY = DAYY00 + DA0Y0Y + DA0Y0Y + DA00YY
           TERMYZ = DAYZ00 + DA0Y0Z + DA0Z0Y + DA00YZ
           TERMZZ = DAZZ00 + DA0Z0Z + DA0Z0Z + DA00ZZ
           DO 400 IREP = 0, MAXLOP
            CHIB = PT(IBTAND(ISYMOP,IREP))
            CSBX = CHIB*SBX
            CSBY = CHIB*SBY
            CSBZ = CHIB*SBZ
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
C   A-A nuclear-attraction Hessian elements
C
            IF (DOAX)
     *         HESSNA(IAX,IAX) = HESSNA(IAX,IAX) + TERMXX
            IF (DOAX.AND.DOAY)
     *         HESSNA(IAX,IAY) = HESSNA(IAX,IAY) + TERMXY
            IF (DOAX.AND.DOAZ)
     *         HESSNA(IAX,IAZ) = HESSNA(IAX,IAZ) + TERMXZ
            IF (DOAY)
     *         HESSNA(IAY,IAY) = HESSNA(IAY,IAY) + TERMYY
            IF (DOAY.AND.DOAZ)
     *         HESSNA(IAY,IAZ) = HESSNA(IAY,IAZ) + TERMYZ
            IF (DOAZ)
     *         HESSNA(IAZ,IAZ) = HESSNA(IAZ,IAZ) + TERMZZ
C
C   B-B nuclear-attraction Hessian elements
C
            IF (DOBX)
     *         HESSNA(IBX,IBX) = HESSNA(IBX,IBX) + TERMXX
            IF (DOBX.AND.DOBY)
     *         HESSNA(IBX,IBY) = HESSNA(IBX,IBY) + SBX*SBY*TERMXY
            IF (DOBX.AND.DOBZ)
     *         HESSNA(IBX,IBZ) = HESSNA(IBX,IBZ) + SBX*SBZ*TERMXZ
            IF (DOBY)
     *         HESSNA(IBY,IBY) = HESSNA(IBY,IBY) + TERMYY
            IF (DOBY.AND.DOBZ)
     *         HESSNA(IBY,IBZ) = HESSNA(IBY,IBZ) + SBY*SBZ*TERMYZ
            IF (DOBZ.AND.DOBZ)
     *         HESSNA(IBZ,IBZ) = HESSNA(IBZ,IBZ) + TERMZZ
C
C   A-B nuclear-attraction Hessian elements
C
            IF (DOAX.AND.DOBX)
     *         HESSNA(IAX,IBX) = HESSNA(IAX,IBX) - FAB*CSBX*TERMXX
            IF (DOAX.AND.DOBY)
     *         HESSNA(IAX,IBY) = HESSNA(IAX,IBY) - CSBY*TERMXY
            IF (DOAX.AND.DOBZ)
     *         HESSNA(IAX,IBZ) = HESSNA(IAX,IBZ) - CSBZ*TERMXZ
            IF (DOAY.AND.DOBX)
     *         HESSNA(IAY,IBX) = HESSNA(IAY,IBX) - CSBX*TERMXY
            IF (DOAY.AND.DOBY)
     *         HESSNA(IAY,IBY) = HESSNA(IAY,IBY) - FAB*CSBY*TERMYY
            IF (DOAY.AND.DOBZ)
     *         HESSNA(IAY,IBZ) = HESSNA(IAY,IBZ) - CSBZ*TERMYZ
            IF (DOAZ.AND.DOBX)
     *         HESSNA(IAZ,IBX) = HESSNA(IAZ,IBX) - CSBX*TERMXZ
            IF (DOAZ.AND.DOBY)
     *         HESSNA(IAZ,IBY) = HESSNA(IAZ,IBY) - CSBY*TERMYZ
            IF (DOAZ.AND.DOBZ)
     *         HESSNA(IAZ,IBZ) = HESSNA(IAZ,IBZ) - FAB*CSBZ*TERMZZ
  400      CONTINUE
          END IF
C
C  two-center case: C = B *****
C
         ELSE IF (ICENTC .EQ. ICENTB) THEN
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
C   A nuclear-attraction gradient elements
C
          IF (DOAX) GRADNA(IAX) = GRADNA(IAX) + DA0X00
          IF (DOAY) GRADNA(IAY) = GRADNA(IAY) + DA0Y00
          IF (DOAZ) GRADNA(IAZ) = GRADNA(IAZ) + DA0Z00
C
C   B nuclear-attraction gradient elements
C
          IF (DOBX) GRADNA(IBX) = GRADNA(IBX) - SBX*DA0X00
          IF (DOBY) GRADNA(IBY) = GRADNA(IBY) - SBY*DA0Y00
          IF (DOBZ) GRADNA(IBZ) = GRADNA(IBZ) - SBZ*DA0Z00
C
C   Second Derivatives:
C
          IF (SECOND) THEN
           DO 500 IREP = 0, MAXLOP
            CHIB = PT(IBTAND(ISYMOP,IREP))
            CSBX = CHIB*SBX
            CSBY = CHIB*SBY
            CSBZ = CHIB*SBZ
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
C    A-A nuclear-attraction Hessian elements
C
            IF (DOAX)
     *         HESSNA(IAX,IAX) = HESSNA(IAX,IAX) + DAXX00
            IF (DOAX.AND.DOAY)
     *         HESSNA(IAX,IAY) = HESSNA(IAX,IAY) + DAXY00
            IF (DOAX.AND.DOAZ)
     *         HESSNA(IAX,IAZ) = HESSNA(IAX,IAZ) + DAXZ00
            IF (DOAY)
     *         HESSNA(IAY,IAY) = HESSNA(IAY,IAY) + DAYY00
            IF (DOAY.AND.DOAZ)
     *         HESSNA(IAY,IAZ) = HESSNA(IAY,IAZ) + DAYZ00
            IF (DOAZ)
     *         HESSNA(IAZ,IAZ) = HESSNA(IAZ,IAZ) + DAZZ00
C
C   B-B nuclear-attraction Hessian elements
C
            IF (DOBX)
     *         HESSNA(IBX,IBX) = HESSNA(IBX,IBX) + DAXX00
            IF (DOBX.AND.DOBY)
     *         HESSNA(IBX,IBY) = HESSNA(IBX,IBY) + SBX*SBY*DAXY00
            IF (DOBX.AND.DOBZ)
     *         HESSNA(IBX,IBZ) = HESSNA(IBX,IBZ) + SBX*SBZ*DAXZ00
            IF (DOBY)
     *         HESSNA(IBY,IBY) = HESSNA(IBY,IBY) + DAYY00
            IF (DOBY.AND.DOBZ)
     *         HESSNA(IBY,IBZ) = HESSNA(IBY,IBZ) + SBY*SBZ*DAYZ00
            IF (DOBZ)
     *         HESSNA(IBZ,IBZ) = HESSNA(IBZ,IBZ) + DAZZ00
C
C   A-B nuclear-attraction Hessian elements
C
            IF (DOAX.AND.DOBX)
     *         HESSNA(IAX,IBX) = HESSNA(IAX,IBX) - FAB*CSBX*DAXX00
            IF (DOAX.AND.DOBY)
     *         HESSNA(IAX,IBY) = HESSNA(IAX,IBY) - CSBY*DAXY00
            IF (DOAX.AND.DOBZ)
     *         HESSNA(IAX,IBZ) = HESSNA(IAX,IBZ) - CSBZ*DAXZ00
            IF (DOAY.AND.DOBX)
     *         HESSNA(IAY,IBX) = HESSNA(IAY,IBX) - CSBX*DAXY00
            IF (DOAY.AND.DOBY)
     *         HESSNA(IAY,IBY) = HESSNA(IAY,IBY) - FAB*CSBY*DAYY00
            IF (DOAY.AND.DOBZ)
     *         HESSNA(IAY,IBZ) = HESSNA(IAY,IBZ) - CSBZ*DAYZ00
            IF (DOAZ.AND.DOBX)
     *         HESSNA(IAZ,IBX) = HESSNA(IAZ,IBX) - CSBX*DAXZ00
            IF (DOAZ.AND.DOBY)
     *         HESSNA(IAZ,IBY) = HESSNA(IAZ,IBY) - CSBY*DAYZ00
            IF (DOAZ.AND.DOBZ)
     *         HESSNA(IAZ,IBZ) = HESSNA(IAZ,IBZ) - FAB*CSBZ*DAZZ00
  500      CONTINUE
          END IF
C
C   three-center case ***** <A|V(C)|B>
C
         ELSE
          IAX  = IPTCNT(NAX,0)
          IAY  = IPTCNT(NAY,0)
          IAZ  = IPTCNT(NAZ,0)
          IBX  = IPTCNT(NBX,0)
          IBY  = IPTCNT(NBY,0)
          IBZ  = IPTCNT(NBZ,0)
          ICX  = IPTCNT(NCX,0)
          ICY  = IPTCNT(NCY,0)
          ICZ  = IPTCNT(NCZ,0)
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
C  A nuclear-attraction gradient elements:
C
          IF (DOAX) GRADNA(IAX) = GRADNA(IAX) + DA0X00
          IF (DOAY) GRADNA(IAY) = GRADNA(IAY) + DA0Y00
          IF (DOAZ) GRADNA(IAZ) = GRADNA(IAZ) + DA0Z00
C
C  B nuclear-attraction gradient elements:
C
          IF (DOBX) GRADNA(IBX) = GRADNA(IBX)
     *                          - SBX*(DA0X00 + DA000X)
          IF (DOBY) GRADNA(IBY) = GRADNA(IBY)
     *                          - SBY*(DA0Y00 + DA000Y)
          IF (DOBZ) GRADNA(IBZ) = GRADNA(IBZ)
     *                          - SBZ*(DA0Z00 + DA000Z)
C
C  C nuclear-attraction gradient elements:
C
          IF (DOCX) GRADNA(ICX) = GRADNA(ICX) + SCX*DA000X
          IF (DOCY) GRADNA(ICY) = GRADNA(ICY) + SCY*DA000Y
          IF (DOCZ) GRADNA(ICZ) = GRADNA(ICZ) + SCZ*DA000Z
C
C   Second Derivatives:
C
            IF (SECOND) THEN
               FAC = D1
               FBC = D1
               IF (NCENTA .EQ. KCENTC) FAC = D2
               IF (NCENTB .EQ. KCENTC) FBC = D2
               DO 600 IREP = 0, MAXLOP
                  CHIB = PT(IBTAND(ISYMOP,IREP))
                  CHIC = PT(IBTAND(ISYMPC,IREP))
                  CSBX = CHIB*SBX
                  CSBY = CHIB*SBY
                  CSBZ = CHIB*SBZ
                  CSCX = CHIC*SCX
                  CSCY = CHIC*SCY
                  CSCZ = CHIC*SCZ
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
C   A-A nuclear-attraction Hessian elements:
C
                  IF (DOAX)
     *               HESSNA(IAX,IAX) = HESSNA(IAX,IAX) + DAXX00
                  IF (DOAX.AND.DOAY)
     *               HESSNA(IAX,IAY) = HESSNA(IAX,IAY) + DAXY00
                  IF (DOAX.AND.DOAZ)
     *               HESSNA(IAX,IAZ) = HESSNA(IAX,IAZ) + DAXZ00
                  IF (DOAY)
     *               HESSNA(IAY,IAY) = HESSNA(IAY,IAY) + DAYY00
                  IF (DOAY.AND.DOAZ)
     *               HESSNA(IAY,IAZ) = HESSNA(IAY,IAZ) + DAYZ00
                  IF (DOAZ)
     *               HESSNA(IAZ,IAZ) = HESSNA(IAZ,IAZ) + DAZZ00
C
C   A-B nuclear-attraction Hessian elements:
C
                  IF (DOAX.AND.DOBX)
     *               HESSNA(IAX,IBX) = HESSNA(IAX,IBX)
     *                               - FAB*CSBX*(DAXX00 + DA0X0X)
                  IF (DOAX.AND.DOBY)
     *               HESSNA(IAX,IBY) = HESSNA(IAX,IBY)
     *                               - CSBY*(DAXY00 + DA0X0Y)
                  IF (DOAX.AND.DOBZ)
     *               HESSNA(IAX,IBZ) = HESSNA(IAX,IBZ)
     *                               - CSBZ*(DAXZ00 + DA0X0Z)
                  IF (DOAY.AND.DOBX)
     *               HESSNA(IAY,IBX) = HESSNA(IAY,IBX)
     *                               - CSBX*(DAXY00 + DA0Y0X)
                  IF (DOAY.AND.DOBY)
     *               HESSNA(IAY,IBY) = HESSNA(IAY,IBY)
     *                               - FAB*CSBY*(DAYY00 + DA0Y0Y)
                  IF (DOAY.AND.DOBZ)
     *               HESSNA(IAY,IBZ) = HESSNA(IAY,IBZ)
     *                               - CSBZ*(DAYZ00 + DA0Y0Z)
                  IF (DOAZ.AND.DOBX)
     *               HESSNA(IAZ,IBX) = HESSNA(IAZ,IBX)
     *                               - CSBX*(DAXZ00 + DA0Z0X)
                  IF (DOAZ.AND.DOBY)
     *               HESSNA(IAZ,IBY) = HESSNA(IAZ,IBY)
     *                               - CSBY*(DAYZ00 + DA0Z0Y)
                  IF (DOAZ.AND.DOBZ)
     *               HESSNA(IAZ,IBZ) = HESSNA(IAZ,IBZ)
     *                               - FAB*CSBZ*(DAZZ00 + DA0Z0Z)
C
C  A-C nuclear-attraction Hessian elements:
C
                  IF (DOAX.AND.DOCX)
     *               HESSNA(IAX,ICX) = HESSNA(IAX,ICX) + FAC*CSCX*DA0X0X
                  IF (DOAX.AND.DOCY)
     *               HESSNA(IAX,ICY) = HESSNA(IAX,ICY) + CSCY*DA0X0Y
                  IF (DOAX.AND.DOCZ)
     *               HESSNA(IAX,ICZ) = HESSNA(IAX,ICZ) + CSCZ*DA0X0Z
                  IF (DOAY.AND.DOCX)
     *               HESSNA(IAY,ICX) = HESSNA(IAY,ICX) + CSCX*DA0Y0X
                  IF (DOAY.AND.DOCY)
     *               HESSNA(IAY,ICY) = HESSNA(IAY,ICY) + FAC*CSCY*DA0Y0Y
                  IF (DOAY.AND.DOCZ)
     *               HESSNA(IAY,ICZ) = HESSNA(IAY,ICZ) + CSCZ*DA0Y0Z
                  IF (DOAZ.AND.DOCX)
     *               HESSNA(IAZ,ICX) = HESSNA(IAZ,ICX) + CSCX*DA0Z0X
                  IF (DOAZ.AND.DOCY)
     *               HESSNA(IAZ,ICY) = HESSNA(IAZ,ICY) + CSCY*DA0Z0Y
                  IF (DOAZ.AND.DOCZ)
     *               HESSNA(IAZ,ICZ) = HESSNA(IAZ,ICZ) + FAC*CSCZ*DA0Z0Z
C
C    B-B nuclear-attraction Hessian elements:
C
                  IF (DOBX)
     *               HESSNA(IBX,IBX) = HESSNA(IBX,IBX) +
     *               (DAXX00 + DA00XX + DA0X0X + DA0X0X)
                  IF (DOBX.AND.DOBY)
     *               HESSNA(IBX,IBY) = HESSNA(IBX,IBY) + SBX*SBY*
     *               (DAXY00 + DA00XY + DA0X0Y + DA0Y0X)
                  IF (DOBX.AND.DOBZ)
     *               HESSNA(IBX,IBZ) = HESSNA(IBX,IBZ) + SBX*SBZ*
     *               (DAXZ00 + DA00XZ + DA0X0Z + DA0Z0X)
                  IF (DOBY)
     *               HESSNA(IBY,IBY) = HESSNA(IBY,IBY) +
     *               (DAYY00 + DA00YY + DA0Y0Y + DA0Y0Y)
                  IF (DOBY.AND.DOBZ)
     *               HESSNA(IBY,IBZ) = HESSNA(IBY,IBZ) + SBY*SBZ*
     *               (DAYZ00 + DA00YZ + DA0Y0Z + DA0Z0Y)
                  IF (DOBZ)
     *               HESSNA(IBZ,IBZ) = HESSNA(IBZ,IBZ) +
     *               (DAZZ00 + DA00ZZ + DA0Z0Z + DA0Z0Z)
C
C   B-C nuclear-attraction Hessian elements:
C
                  IF (DOBX.AND.DOCX) HESSNA(IBX,ICX) = HESSNA(IBX,ICX)
     *                  - FBC* CSBX*CSCX*(DA0X0X + DA00XX)
                  IF (DOBX.AND.DOCY) HESSNA(IBX,ICY) = HESSNA(IBX,ICY)
     *                  - CSBX*CSCY*(DA0X0Y + DA00XY)
                  IF (DOBX.AND.DOCZ) HESSNA(IBX,ICZ) = HESSNA(IBX,ICZ)
     *                  - CSBX*CSCZ*(DA0X0Z + DA00XZ)
                  IF (DOBY.AND.DOCX) HESSNA(IBY,ICX) = HESSNA(IBY,ICX)
     *                  - CSBY*CSCX*(DA0Y0X + DA00XY)
                  IF (DOBY.AND.DOCY) HESSNA(IBY,ICY) = HESSNA(IBY,ICY)
     *                  - FBC* CSBY*CSCY*(DA0Y0Y + DA00YY)
                  IF (DOBY.AND.DOCZ) HESSNA(IBY,ICZ) = HESSNA(IBY,ICZ)
     *                  - CSBY*CSCZ*(DA0Y0Z + DA00YZ)
                  IF (DOBZ.AND.DOCX) HESSNA(IBZ,ICX) = HESSNA(IBZ,ICX)
     *                  - CSBZ*CSCX*(DA0Z0X + DA00XZ)
                  IF (DOBZ.AND.DOCY) HESSNA(IBZ,ICY) = HESSNA(IBZ,ICY)
     *                  - CSBZ*CSCY*(DA0Z0Y + DA00YZ)
                  IF (DOBZ.AND.DOCZ) HESSNA(IBZ,ICZ) = HESSNA(IBZ,ICZ)
     *                  - FBC*CSBZ*CSCZ*(DA0Z0Z + DA00ZZ)
C
C   C-C nuclear-attraction Hessian elements:
C
                  IF (DOCX)
     *               HESSNA(ICX,ICX) = HESSNA(ICX,ICX) + DA00XX
                  IF (DOCX.AND.DOCY)
     *               HESSNA(ICX,ICY) = HESSNA(ICX,ICY) + SCX*SCY*DA00XY
                  IF (DOCX.AND.DOCZ)
     *               HESSNA(ICX,ICZ) = HESSNA(ICX,ICZ) + SCX*SCZ*DA00XZ
                  IF (DOCY)
     *               HESSNA(ICY,ICY) = HESSNA(ICY,ICY) + DA00YY
                  IF (DOCY.AND.DOCZ)
     *               HESSNA(ICY,ICZ) = HESSNA(ICY,ICZ) + SCY*SCZ*DA00YZ
                  IF (DOCZ)
     *               HESSNA(ICZ,ICZ) = HESSNA(ICZ,ICZ) + DA00ZZ
  600          CONTINUE
            END IF
         END IF
  250 CONTINUE
      RETURN
      END
