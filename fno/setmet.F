      
      SUBROUTINE SETMET(UHF)
      IMPLICIT NONE
#include "flags.h"
C     COMMON BLOCKS
      INTEGER IFLAGS(100)
      COMMON/FLAGS/IFLAGS
      INTEGER IFLAGS2(500)
      COMMON/FLAGS2/IFLAGS2
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer naobasis,naobas(8)
      common/aoinfo/naobasis,naobas
      logical rohf,qrhf,semi,natorb,anaderiv
      common/fnoparam/rohf,qrhf,semi,natorb,anaderiv
C     INPUT/OUTPUT VARIABLES
      INTEGER UHF
C------------------------------------------------------------------------
      
 3000 FORMAT('  The reference state is a ROHF wave function.')
 3001 FORMAT('  The reference state uses semi-canonical orbitals.')
 3004 FORMAT('  Untruncated virtual natural orbitals will be formed.')
 3005 FORMAT('  The reference state is a QRHF wave function.')
 3006 FORMAT('  The reference state is a general non-HF wave function.')
 3007 FORMAT('  WARNING!!! Symmetries of the reference may not be',
     &       ' preserved!')
C

      ROHF=.FALSE.
      SEMI=.FALSE.
      QRHF=.FALSE.
      ANADERIV=(IFLAGS(h_IFLAGS_deriv_lev).GT.0).OR.
     &   (IFLAGS(h_IFLAGS_props).GT.0)

C
C CHECK IF WE ARE DEALING WITH A NON-HF CASE 
C
      NATORB=.FALSE.
      IF (IFLAGS2(h_IFLAGS2_natural).EQ.1) then
         NATORB=.TRUE.
         write(6,3004)
      endif

      IF ((IFLAGS(h_IFLAGS_reference).eq.2).or.
     &   (IFLAGS(h_IFLAGS_qrhf_gen).ne.0)) then

C     ROHF AND QRHF cases
         write(6,3006)
         IF(IFLAGS(h_IFLAGS_reference).EQ.2) THEN
            ROHF=.TRUE.
            WRITE(6,3000)
            UHF=1
         ENDIF

         IF (IFLAGS(h_IFLAGS_qrhf_gen).ne.0) then
            qrhf=.true.
            write(6,3005)
            write(6,*) 'QRHF FNOs not currently supported.'
            call errex
            write(6,3007)
            UHF=1
         ENDIF
      ELSE IF (IFLAGS(h_IFLAGS_nonhf).eq.1) then
C     Other NONHF cases
         write(6,3006)
         write(6,3007)
      ENDIF
      
      IF (IFLAGS(h_IFLAGS_orbitals).eq.1) THEN
         SEMI=.TRUE.
         WRITE(6,3001)
      ENDIF
      
      CALL GETREC(20,'JOBARC','NBASTOT',1,NAOBASIS)
      CALL GETREC(20,'JOBARC','NUMBASIR',NIRREP,NAOBAS)

      RETURN
      END
