      SUBROUTINE QU0PRI(QUADM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "mxcent.par"
#include "baslims.par"
      CHARACTER*1 CHRXYZ(-3:3)
      PARAMETER (DEBYE = 2.541 770 D00)
      LOGICAL SKIP, NODC, NODV, TEST
      COMMON /CBIDIP/ IPRINT, SKIP, NODC, NODV, TEST
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
      DIMENSION QUADM(6)
C
      DATA CHRXYZ /'z','y','x',' ','X','Y','Z'/
C
      WRITE (LUPRI,'(31X,A/)') 'au            '
      IND=0
      DO 100 I = 1, 3
       DO 100 J=I,3
        IND=IND+1
        IF (ISYTYP(I).EQ.ISYTYP(J)) THEN
         WRITE (LUPRI,'(20X,A,A,F16.8)') CHRXYZ(-I),CHRXYZ(-J),
     *             QUADM(IND)
        END IF
100    CONTINUE
c      WRITE (LUPRI,'(/A,F12.8,A)')
c     *   ' Conversion factor used: 1 Debye = ',DEBYE,' a.u. '
      WRITE (LUPRI,'(/)')
      RETURN
      END
