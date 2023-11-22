#include <flags.h>
      SUBROUTINE DRIVE_RPA(ICORE,MAXCOR,IUHF)
C
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR)
      LOGICAL RCCD_EOM
C
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /EOM_METHOD/RCCD_EOM
C
#include "flags.com"
#include "machsp.com"

C If any of the r/dr-CCD EOM flags are set and RPA is also run, then
C write the RPA T2 vectors to lists and set a flag to inform the r/dr-CCD
C codes that T2 is already available.

      RCCD_EOM = ((IFLAGS(h_IFLAGS_excite) .EQ. 13) .OR.
     +            (IFLAGS(h_IFLAGS_excite) .EQ. 14) .OR.
     +            (IFLAGS(h_IFLAGS_excite) .EQ. 15) .OR.
     +            (IFLAGS(h_IFLAGS_excite) .EQ. 16) .OR.
     +            (IFLAGS(h_IFLAGS_excite) .EQ. 17) .OR.
     +            (IFLAGS(h_IFLAGS_excite) .EQ. 18) .OR.
     +            (IFLAGS(h_IFLAGS_excite) .EQ. 19) .OR.
     +            (IFLAGS(h_IFLAGS_excite) .EQ. 20))

      IF (IUHF.EQ.0)THEN
C
          NSIZ=IRPDPD(1,ISYTYP(1,25))
C
          I000=1
          I010=I000+2*IINTFP*NSIZ
          MXCOR=MAXCOR-I010+1
          CALL DO_RHF_RPA(ICORE(I010),MXCOR,ICORE(I000),NSIZ,IUHF)
C
      ELSEIF (IUHF.EQ.1) THEN

          NSIZ=IRPDPD(1,ISYTYP(1,19))+IRPDPD(1,ISYTYP(1,20))

          I000=1
          I010=I000+2*IINTFP*NSIZ
          MXCOR=MAXCOR-I010+1
          CALL DO_UHF_RPA(ICORE(I010),MXCOR,ICORE(I000),NSIZ,IUHF)

      ENDIF

      IF (RCCD_EOM) CALL PROCESS_RPAT2(ICORE,MAXCOR,IUHF)
C
CSS      I000 = 1
CSS      call modf_amps(ICORE(I000),MAXCOR,IUHF,0,.FALSE.,'T',"SAVE",
CSS     &               "AMPS")
C
C Change the refrence to RHF to compute the T1 energy contribution 
C corresponding to direct RPA.
C
      iflags(h_IFLAGS_reference) = 1
      Iuhf = 0
CSSS      call Do_t1crs_2RPA(ICORE(I000), MAXCOR, IUHF)
C
      RETURN
      END
