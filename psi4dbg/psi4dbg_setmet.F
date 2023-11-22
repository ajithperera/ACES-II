#include "flags.h"
      SUBROUTINE PSI4DBG_SETMET()

      IMPLICIT INTEGER(A-Z)

      LOGICAL NONHF
      LOGICAL CIS,EOM
      LOGICAL PCCD,CCD,LCCD,MBPT2

#include "jodaflags.com"
#include "machsp.com"

      COMMON/METH/PCCD,CCD,LCCD,MBPT2
      COMMON/EXCITE/EOM,CIS
      COMMON/NHFREF/NONHF
      COMMON/EXTRAPO/RLE 
C
      PCCD    = .FALSE.
      CCD     = .FALSE.
      LCCD    = .FALSE.
      EOM     = .FALSE.
      CIS     = .FALSE.
      NONHF   = .FALSE.
      RLE     = .FALSE.
      MBPT2   = .FALSE.

      IF(IFLAGS(h_IFLAGS_nonhf)  .EQ.1) NONHF = .TRUE.
      IF(IFLAGS(h_IFLAGS_excite) .EQ.3) EOM   = .TRUE.
      IF(IFLAGS(h_IFLAGS_excite) .EQ.5) CIS   = .TRUE.

      IF(IFLAGS(h_IFLAGS_cc_extrapol) .Ne.3) RLE  = .TRUE.

      IF (IFLAGS(h_IFLAGS_calclevel)  .EQ.52) PCCD   = .TRUE.
      IF (IFLAGS(h_IFLAGS_calclevel)  .EQ.53) LCCD   = .TRUE.
      IF (IFLAGS(h_IFLAGS_calclevel)  .EQ.54)  CCD   = .TRUE.
      IF (IFLAGS(h_IFLAGS_calclevel)  .EQ.1) MBPT2   = .TRUE.

C Here LCCD and CCD stands for oo-LCCD and oo-CCD respectively.

      IF (PCCD) THEN
          Write(6,*) 
          Write(6,"(a,a)") " The left and right Pair coupled cluster",
     +                     " doubles (pCCD) equations are solved."
          Write(6,*)
      ELSE IF (LCCD) THEN
          Write(6,*) 
          Write(6,"(a,a)") " The left and right linear coupled cluster",
     +                     " doubles (LCCD) equations are solved."
          Write(6,*)
      ELSE IF (CCD) THEN
          Write(6,*) 
          Write(6,"(a,a)") " The left and right coupled cluster",
     +                     " doubles (CCD) equations are solved."
          Write(6,*) 
      ELSE IF (MBPT2) THEN
          Write(6,*) 
          Write(6,"(a,a)") " The left and right MBPT(2) equations",
     +                     " are solved."
          Write(6,*) 
      ELSE
          Write(6,*) 
          Write(6,"(2a)") " Incorrect choice of method. Choose pCCD",
     +                    " OO-CCD or OO-LCCD instead"
          Call Errex
      ENDIF 

      RETURN
      END
