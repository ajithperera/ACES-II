      SUBROUTINE DP0PRI(DIPOLE)
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
      DIMENSION DIPOLE(3)
C
      DATA CHRXYZ /'z','y','x',' ','X','Y','Z'/
C
      IF (NAXREP(0) .EQ. 0) THEN
         WRITE (LUPRI,'(20X,A)') ' All components zero by symmetry'
      ELSE
         WRITE (LUPRI,'(31X,A/)') 'au             Debye'
         DO 100 I = 1, 3
            IF (ISYTYP(I) .EQ. 0) THEN
               WRITE (LUPRI,'(20X,A,2F16.8)') CHRXYZ(-I),
     *              DIPOLE(I), DEBYE*DIPOLE(I)
            END IF
  100    CONTINUE
      END IF
      WRITE (LUPRI,'(/A,F12.8,A)')
     *   ' Conversion factor used: 1 Debye = ',DEBYE,' a.u. '
      WRITE (LUPRI,'(/)')
      RETURN
      END
