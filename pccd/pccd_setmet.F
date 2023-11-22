#include "flags.h"
      SUBROUTINE PCCD_SETMET()

      IMPLICIT INTEGER(A-Z)

      LOGICAL NONHF_REF,SYMMETRY
      LOGICAL CIS,EOM
      LOGICAL PCCD,CCD,LCCD
      LOGICAL Opt_orbs
      LOGICAL ONEP_ONLY 
      LOGICAL DIIS 

#include "flags.com"
#include "flags2.com"
#include "machsp.com"

      COMMON/CALC/PCCD,CCD,LCCD
      COMMON/EXCITE/EOM,CIS
      COMMON/NHFREF/NONHF_REF
      COMMON/EXTRAPO/DIIS 
      COMMON/SYMM/SYMMETRY
      COMMON/ORBR_HESS/ONEP_ONLY
C
      PCCD     = .FALSE.
      CCD      = .FALSE.
      LCCD     = .FALSE.
      EOM      = .FALSE.
      CIS      = .FALSE.
      NONHF_REF= .FALSE.
      RLE      = .FALSE.
      SYMMETRY = .FALSE. 
      DIIS     = .FALSE. 

C CCD models are not compatible with NONHF refernce key-word simply because NONHF
C wants to inculde singles indirectly. 

      IF(IFLAGS(h_IFLAGS_nonhf)   .EQ.1) NONHF_REF = .TRUE.
      IF(IFLAGS(h_IFLAGS_symmetry).EQ.2) SYMMETRY  = .TRUE.

      IF(IFLAGS(h_IFLAGS_excite) .EQ.3) EOM   = .TRUE.
      IF(IFLAGS(h_IFLAGS_excite) .EQ.5) CIS   = .TRUE.

      IF(IFLAGS2(h_IFLAGS2_oo_diis) .Ne.0) DIIS = .TRUE.

      IF (IFLAGS(h_IFLAGS_calclevel)  .EQ.52) PCCD   = .TRUE.
      IF (IFLAGS(h_IFLAGS_calclevel)  .EQ.53) LCCD   = .TRUE.
      IF (IFLAGS(h_IFLAGS_calclevel)  .EQ.54)  CCD   = .TRUE.

      ONEP_ONLY = .TRUE. 

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
      ELSE
          Write(6,*) 
          Write(6,"(2a)") " Incorrect choice of method. Choose pCCD",
     +                    " OO-CCD or OO-LCCD instead"
          Call Errex
      ENDIF 

      RETURN
      END
