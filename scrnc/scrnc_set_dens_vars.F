#include "flags.h"
      Subroutine Scrnc_set_dens_vars(IUHF)
    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER*80 NAME,NAME1,NAME2,NAME3,NAME4,NAME5,NAME6,NAME7

      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD_DUPL,UCC
      LOGICAL DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,
     &             TRIP2
      LOGICAL EOM,RELAXED,CIS,PROPS
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD_DUPL,
     &            UCC
      COMMON/DERIV/DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,
     &             TRIP2

#include "jodaflags.com"
C
C These variables are set primarily for legacy and the current code
C is applicable only to CCSD,RCCD,DRCCD response densities.
  
      MBPT2       = .FALSE.
      MBPT3       = .FALSE.
      M4DQ        = .FALSE.
      M4SDQ       = .FALSE.
      M4SDTQ      = .FALSE.
      CCD         = .FALSE.
      QCISD       = .FALSE.
      CCSD_DUPL   = .FALSE.
      UCC         = .FALSE.
      NONHF       = .FALSE.
      ROHF        = .FALSE.
      SEMI        = .FALSE.
      CANON       = .FALSE.
      TRIP1       = .FALSE.
      TRIP2       = .FALSE.

      IF(IFLAGS(h_IFLAGS_calc).EQ.0) THEN
         IF(IFLAGS(h_IFLAGS_excite).EQ.1)THEN
            CIS=.TRUE.
            NAME2='TDA'
            write(6,3002)NAME2
         ENDIF
      ELSEIF (IFLAGS2(h_IFLAGS2_program) .EQ. 3) THEN
         CCSD_DUPL = .TRUE.
         EOM=.TRUE.
         write(6,3021)
      ELSEIF (IFLAGS(h_IFLAGS_excite).EQ.17
     +        .OR.IFLAGS(h_IFLAGS_excite).EQ.18)  THEN
         CCSD_DUPL=.TRUE.
         NAME3='CCSD'
         write(6,3020)
         EOM=.TRUE.
      ELSEIF(IFLAGS(h_IFLAGS_calc).EQ.1) THEN
         MBPT2=.TRUE.
         NAME='MBPT(2)'
         write(6,3000)NAME
      ELSE IF(IFLAGS(h_IFLAGS_calc).EQ.2) THEN
         MBPT3=.TRUE.
         NAME='MBPT(3)'
         write(6,3000)NAME
      ELSE IF(IFLAGS(h_IFLAGS_calc).EQ.3) THEN
         M4SDQ=.TRUE.
         NAME1='SDQ-MBPT(4)'
         write(6,3001) NAME1
      ELSE IF(IFLAGS(h_IFLAGS_calc).EQ.4) THEN
         M4SDTQ=.TRUE.
         TRIP1=.TRUE.
         NAME='MBPT(4)'
         write(6,3000)NAME
      ELSE IF(IFLAGS(h_IFLAGS_calc).EQ.7) THEN
         UCC=.TRUE.
         NAME5='UCCSD(4)'
         write(6,3010)NAME5
      ELSE IF(IFLAGS(h_IFLAGS_calc).EQ.8) THEN
         CCD=.TRUE.
         NAME2='CCD'
         write(6,3002)NAME2
      ELSE IF(IFLAGS(h_IFLAGS_calc).EQ.9) THEN
         UCC=.TRUE.
         TRIP1=.TRUE.
         NAME7='UCC(4)'
         write(6,3012)NAME7
      ELSE IF(IFLAGS(h_IFLAGS_calc).EQ.10) THEN
         CCSD_DUPL=.TRUE.
         NAME3='CCSD'
         write(6,3003)NAME3
      ELSE IF(IFLAGS(h_IFLAGS_calc).EQ.11) THEN
         CCSD_DUPL=.TRUE.
         NAME6='CCSD+T(CCSDT)'
         TRIP1=.TRUE.
         write(6,3011) NAME6
      ELSE IF(IFLAGS(h_IFLAGS_calc).EQ.21) THEN
         QCISD=.TRUE.
         NAME5='QCISD(T)'
         TRIP1=.TRUE.
         TRIP2=.TRUE.
         write(6,3010) NAME5
      ELSE IF(IFLAGS(h_IFLAGS_calc).EQ.22) THEN
         CCSD_DUPL=.TRUE.
         NAME='CCSD(T)'
         TRIP1=.TRUE.
         TRIP2=.TRUE.
         write(6,3000) NAME
      ELSE IF(IFLAGS(h_IFLAGS_calc).EQ.23 .OR. 
     +        IFLAGS(h_IFLAGS_calc).EQ. 25) THEN
         QCISD=.TRUE.
         NAME4='QCISD'
         write(6,3004)NAME4
         Write(6,"(a,a)") "  WARNING! CISD or QCISD density matrices",
     +                    " are not fully tested!!!"
      ELSE IF(IFLAGS(h_IFLAGS_calc).EQ.42) THEN
         CCSD_DUPL=.TRUE.
         NAME5='aCCSD(T)'
         TRIP1=.TRUE.
         TRIP2=.TRUE.
         LTRP=.TRUE.
         write(6,3010) NAME5
      ENDIF

      DENS=.TRUE.
      RELAXED = (IFLAGS(h_IFLAGS_density).EQ.0)

      GRAD  = (IFLAGS2(h_IFLAGS2_grad_calc).EQ.1)
      PROPS = (IFLAGS(h_IFLAGS_props)      .EQ.1)

      IF(IFLAGS(h_IFLAGS_qrhf_gen).NE.0)THEN
         QRHF =.TRUE.
         NONHF=.TRUE.
      ENDIF 

      IF(IFLAGS(h_IFLAGS_ref).EQ.2) THEN
C    
         ROHF=.TRUE.
         WRITE(6,3007)
         IUHF=1
         IF(IFLAGS(h_IFLAGS_orbitals).EQ.1) THEN
            SEMI=.TRUE.
            WRITE(6,3008)
         ENDIF
      ENDIF
C
      IF(IFLAGS(h_IFLAGS_pert_orb).EQ.1) THEN
         CANON=.TRUE.
         WRITE(6,3009)
      ENDIF

 3000 FORMAT('  ',A7,' density and intermediates are calculated.')
 3001 FORMAT('  ',A11,' density and intermediates are calculated.')
 3002 FORMAT('  ',A,' density and intermediates are calculated.')
 3003 FORMAT('  ',A4,' density and intermediates are calculated.')
 3004 FORMAT('  ',A5,' density and intermediates are calculated.')
 3005 FORMAT('  The reference state is a QRHF wave functions.')
 3006 FORMAT('  The reference state is a non HF wave functions.')
 3007 FORMAT('  The reference state is a ROHF wave function.')
 3008 FORMAT('  Semi-canonical orbitals are used.')
 3009 FORMAT('  The perturbed orbitals are chosen canonical.')
 3010 FORMAT('  ',A8,' density and intermediates are calculated.')
 3011 FORMAT('  ',A12,' density and intermediates are calculated.')
 3012 FORMAT('  ',A6,' density and intermediates are calculated.')
 3013 FORMAT('  ',A4,' unrelaxed density is evaluated.')
 3020 FORMAT(' RCCD or DRCCD density and intermediates are calculated.')
 3021 FORMAT(' ACES3 density and intermediates are calculated.')

      RETURN
      END
