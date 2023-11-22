#include "flags.h"
      SUBROUTINE RCC_SETMET

      IMPLICIT INTEGER(A-Z)

      LOGICAL NONHF
      LOGICAL CIS,EOM
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD,RLE,RPA_T2VECS

#include "jodaflags.com"
#include "machsp.com"

      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
      COMMON/EXCITE/EOM,CIS
      COMMON/NHFREF/NONHF
      COMMON/EXTRAPO/RLE 
      COMMON/T2_SOURCE/RPA_T2VECS
C
      MBPT2   = .FALSE.
      MBPT3   = .FALSE.
      M4DQ    = .FALSE.
      M4SDQ   = .FALSE.
      M4SDTQ  = .FALSE.
      CCD     = .FALSE.
      CC2     = .FALSE.
      QCISD   = .FALSE.
      CCSD    = .FALSE.
      UCC     = .FALSE.
      NONHF   = .FALSE.
      TRIPIT  = .FALSE.
      TRIPNI  = .FALSE.
      TRIPNI1 = .FALSE.
      RCCD    = .FALSE. 
      DRCCD   = .FALSE.
      EOM     = .FALSE.
      CIS     = .FALSE.
      NONHF   = .FALSE.
      RLE     = .FALSE.

      IF(IFLAGS(h_IFLAGS_nonhf)  .EQ.1) NONHF = .TRUE.
      IF(IFLAGS(h_IFLAGS_excite) .EQ.3) EOM   = .TRUE.
      IF(IFLAGS(h_IFLAGS_excite) .EQ.5) CIS   = .TRUE.

      IF(IFLAGS(h_IFLAGS_cc_extrapol) .Ne.3) RLE  = .TRUE.

      IF (IFLAGS(h_IFLAGS_calclevel)  .EQ.48) RCCD   = .TRUE.
      IF (IFLAGS(h_IFLAGS_calclevel)  .EQ.49) DRCCD  = .TRUE.

      CALL GETREC(0,"JOBARC","RPAT2VEC",LENGTH,D)
      IF (LENGTH .GE. 0) THEN
         WRITE(6,*)
         WRITE(6,"(a,a)") " The RPA module has been run and RPA T2",
     &                    " vectors are read and r/dr-CCD iterations"
         WRITE(6,"(a,a)") " are skipped after constructing the",
     &                    " r/dr-CCD W(mb,ej) intermediate."
         WRITE(6,*)
         RPA_T2VECS = .TRUE.
      ENDIF 

      IF (.NOT. RPA_T2VECS) THEN
      IF (RCCD) THEN
          Write(6,"(a,a)") " Ring coupled cluster doubles (rCCD)",
     +                     " equations are solved."
          Write(6,*)
      ELSE IF (DRCCD) THEN
          Write(6,"(a,a)") " Direct Ring coupled cluster doubles",
     +                     " (drCCSD) equations are solved."
          Write(6,*)
      ENDIF 
      ENDIF 

      RETURN
      END
