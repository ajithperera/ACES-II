#include "flags.h"
      SUBROUTINE RCL_SETMET(HBAR_ONLY)

      IMPLICIT INTEGER(A-Z)

      LOGICAL NONHF
      LOGICAL CIS,EOM,RPA,EOM_SDRCCD,EOM_SRCCD,EOM_SFDRCCD,EOM_SFRCCD
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD,RLE,EOM_SDXRCCD,EOM_SFDXRCCD

      LOGICAL TRIPIT,TRIPNI,TRIPNI1,EOM_SRXCCD,EOM_SFRXCCD,
     &        HBAR_ONLY
     
#include "jodaflags.com"
#include "machsp.com"

      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
      COMMON/EXCITE/EOM,CIS,RPA,EOM_SRCCD,EOM_SDRCCD,EOM_SFRCCD,
     &              EOM_SFDRCCD,EOM_SDXRCCD,EOM_SFDXRCCD
      COMMON/NHFREF/NONHF
      COMMON/EXTRAPO/RLE 
C
      MBPT2       = .FALSE.
      MBPT3       = .FALSE.
      M4DQ        = .FALSE.
      M4SDQ       = .FALSE.
      M4SDTQ      = .FALSE.
      CCD         = .FALSE.
      CC2         = .FALSE.
      QCISD       = .FALSE.
      CCSD        = .FALSE.
      UCC         = .FALSE.
      NONHF       = .FALSE.
      TRIPIT      = .FALSE.
      TRIPNI      = .FALSE.
      TRIPNI1     = .FALSE.
      RCCD        = .FALSE. 
      DRCCD       = .FALSE.
      EOM         = .FALSE.
      CIS         = .FALSE.
      NONHF       = .FALSE.
      RLE         = .FALSE.
      RPA         = .FALSE. 
      EOM_SFDRCCD = .FALSE.
      EOM_SFRCCD  = .FALSE.
      EOM_SDRCCD  = .FALSE.
      EOM_SRCCD   = .FALSE.
      EOM_SRXCCD  = .FALSE.
      EOM_SFRXCCD = .FALSE.

      IF(IFLAGS(h_IFLAGS_nonhf)  .EQ.1)  NONHF        = .TRUE.
      IF(IFLAGS(h_IFLAGS_excite) .EQ.3)  EOM          = .TRUE.
      IF(IFLAGS(h_IFLAGS_excite) .EQ.5 .OR. 
     +   IFLAGS(h_IFLAGS_excite) .EQ.1)  CIS          = .TRUE.
      IF(IFLAGS(h_IFLAGS_excite) .EQ.2)  RPA          = .TRUE.
      IF(IFLAGS(h_IFLAGS_excite) .EQ.17) EOM_SFRCCD   = .TRUE.
      IF(IFLAGS(h_IFLAGS_excite) .EQ.18) EOM_SFDRCCD  = .TRUE.
      IF(IFLAGS(h_IFLAGS_excite) .EQ.13) EOM_SDRCCD   = .TRUE.
      IF(IFLAGS(h_IFLAGS_excite) .EQ.15) EOM_SRCCD    = .TRUE.
      IF(IFLAGS(h_IFLAGS_excite) .EQ.19) EOM_SDXRCCD  = .TRUE.
      IF(IFLAGS(h_IFLAGS_excite) .EQ.20) EOM_SFDXRCCD = .TRUE.

      IF(IFLAGS(h_IFLAGS_cc_extrapol) .Ne.3) RLE  = .TRUE.

      IF (IFLAGS(h_IFLAGS_calclevel)  .EQ.48) RCCD   = .TRUE.
      IF (IFLAGS(h_IFLAGS_calclevel)  .EQ.49) DRCCD  = .TRUE.

      IF (RPA .AND. RCCD)  EOM_SFRCCD  = .TRUE.
      IF (RPA .AND. DRCCD) EOM_SFDRCCD = .TRUE.

      IF (HBAR_ONLY) RETURN
 
      IF (RCCD) THEN
          Write(6,"(a,a)") " Ring coupled cluster doubles (rCCD)",
     +                     " lambda equations are solved."
          Write(6,*)
      ELSE IF (DRCCD) THEN
          Write(6,"(a,a)") " Direct Ring coupled cluster doubles",
     +                     " (drCCSD) lambda equations are solved."
          Write(6,*)
      ENDIF 
      RETURN
      END
