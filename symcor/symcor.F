#include "flags.h"
C
C THIS PROGRAM DETERMINES SYMMETRY ADAPTED COORDINATES WHICH
C ARE USED FOR NUMERICAL DERIVATIVE CALCULATIONS.  IT
C CAN EASILY BE ADAPTED TO CONSTRUCT OTHER SYMMETRY ADAPTED
C QUANTITIES, SUCH AS BASIS FUNCTIONS.  IT CURRENTLY WORKS FOR
C ALL POINT GROUPS OTHER THAN THOSE WHICH CONTAIN COMPLEX
C REPRESENTATIONS (CN, SN, CNH (N>2) AND T AND TH).  FOR THESE
C "DANGEROUS GROUPS", THE ABELIAN SUBGROUP IS USED.
C WRITTEN BY J.F. STANTON, 1991-1992.

c The entire operating procedure of symcor is explained in README.

      SUBROUTINE SYMCOR(ICORE,ICRSIZ)
c
      implicit none
      INTEGER ICRSIZ
      INTEGER ICORE(ICRSIZ)
      LOGICAL FDS_OF_VECTRS , FDS_OF_MATELMS
C
#include "control.com"
#include "parallel_aces.com"
#include "jodaflags.com"

C The symcor_onedim handles numerical derivatives of energy, gradients,
C dipole moments etc. We can loosly say that it handles the derivative
C of vectors and the symcor_muldim handles the deriavtives of matrices.
C Personaly, I think their should not be such a distinction and a one
C general routine can handle both but I am not here only to correct past
C mistakes. Ajith Perera 04/07.
C
      FDS_OF_MATELMS = (Iflags2(h_IFLAGS2_prep_vibron) .EQ. 1
     $     .OR. Iflags2(h_IFLAGS2_scan_normal) .EQ. 1
     $     .OR. Iflags2(h_IFLAGS2_geo_sample) .EQ. 1
     $     .OR. Iflags2(h_IFLAGS2_gen_quartic) .EQ. 1 )

      IF (FDS_OF_MATELMS) THEN
          CALL SYMCOR_MULDIM(ICORE, ICRSIZ)
      ELSE

          CALL SYMCOR_ONEDIM(ICORE, ICRSIZ)
      ENDIF

      RETURN
      END

