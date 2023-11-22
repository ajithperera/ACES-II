      SUBROUTINE CLCNRXLFT(NROWXLFT)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Argument.
C-----------------------------------------------------------------------
      INTEGER NROWXLFT
C-----------------------------------------------------------------------
C     Common block variables.
C-----------------------------------------------------------------------
      DOUBLE PRECISION PC
      INTEGER DSTRT,NTAP,LU2,NRSS,NUCZ,ITAG,MAXLOP,MAXLOT,KMAX,NMAX,
     &        KHKT,MULT,ISYTYP,ITYPE,AND,OR,EOR,NPARSU,NPAR,MULNUC,
     &        NHKT,MUL,NUCO,NRCO,JSTRT,NSTRT,MST,JRS
C-----------------------------------------------------------------------
C     Local variables.
C-----------------------------------------------------------------------
      INTEGER IA,IB,ITOT,NC,ND
C-----------------------------------------------------------------------
C     Parameters.
C-----------------------------------------------------------------------
#include "baslims.par"
#include "mxatms.par"
C-----------------------------------------------------------------------
      COMMON /INDX/ PC(512),DSTRT(8,MXCBF),NTAP,LU2,NRSS,NUCZ,ITAG,
     & MAXLOP,MAXLOT,KMAX,NMAX,KHKT(7),MULT(8),ISYTYP(3),ITYPE(7,28),
     & AND(8,8),OR(8,8),EOR(8,8),NPARSU(8),NPAR(8),MULNUC(Mxatms),
     &  NHKT(MXTNSH),  MUL(MXTNSH),NUCO(MXTNSH),NRCO(MXTNSH),
     & JSTRT(MXTNSH),NSTRT(MXTNSH), MST(MXTNSH), JRS(MXTNSH)
C
C-----------------------------------------------------------------------
CJDW 1/28/97.
C     Block of code to determine maximum leading dimension of XLFT.
C     Previously this was hardcoded to 4095, which is not enough for
C     the extremely common case of 14 d primitives, for example.
C
      NROWXLFT = 0
      DO 100 IA=1,KMAX
      DO  90 IB=1,IA
C
      ITOT=0
      DO  80 NC=1,NHKT(IA)
      DO  70 ND=1,NHKT(IB)
C
      ITOT = ITOT + NC + ND - 1
   70 CONTINUE
   80 CONTINUE
      ITOT = ITOT * NUCO(IA) * NUCO(IB)
C
      NROWXLFT = MAX(NROWXLFT,ITOT)
   90 CONTINUE
  100 CONTINUE
C
c     write(6,*) kmax,nhkt,nuco,itot
c     write(6,*) '  @CLCNRXLFT-I, NROWXLFT ',NROWXLFT
      IF(NROWXLFT .LT. 4095) NROWXLFT = 4095
cYAU      write(6,*) '  @CLCNRXLFT-I, NROWXLFT ',NROWXLFT
C-----------------------------------------------------------------------
      RETURN
      END
