      
      SUBROUTINE SETMET(IUHF)
      IMPLICIT INTEGER(A-Z)
      CHARACTER*20 NAME
      CHARACTER*11 NAME1
      CHARACTER*3 NAME2
      CHARACTER*4 NAME3
      CHARACTER*5 NAME4
      CHARACTER*8 NAME5
      CHARACTER*12 NAME6
      CHARACTER*6 NAME7
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      LOGICAL DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,TRIP2,
     &     GABCD,RELAXED,TRULY_NONHF, PROPS, LTRP
      LOGICAL CIS,EOM
      LOGICAL pCCD,NTORB
      COMMON/EXCITE/CIS,EOM
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON/DERIV/DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,
     &     TRIP2,GABCD,RELAXED,TRULY_NONHF
      COMMON/LTRIP/LTRP
#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */
#include "flags.com"
#include "flags2.com"

      EQUIVALENCE(METHOD,IFLAGS(h_IFLAGS_calc))
C      
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
 3014 FORMAT('  Direct algorithm used for all contributions ',
     &     'involving G(AB,CD).')
 3020 FORMAT('  EOM-',A4,' density and intermediates are calculated.')
 3021 FORMAT('  ACES3 density and intermediates are calculated.')
C     
      MBPT2=.FALSE.
      MBPT3=.FALSE.
      M4DQ=.FALSE.
      M4SDQ=.FALSE.
      M4SDTQ=.FALSE.
      CCD=.FALSE.
      QCISD=.FALSE.
      CCSD=.FALSE.
      UCC=.FALSE.
      NONHF=.FALSE.
      ROHF=.FALSE.
      SEMI=.FALSE.
      CANON=.FALSE.
      TRIP1=.FALSE.
      TRIP2=.FALSE.
      CIS=.FALSE.
      EOM=.FALSE.
      GABCD=.FALSE.
      LTRP=.FALSE.
C     
      IF(METHOD.EQ.0) THEN
         IF(IFLAGS(h_IFLAGS_excite).EQ.1)THEN
            CIS=.TRUE.
            NAME2='TDA'
            write(6,3002)NAME2
         ENDIF
C     SG 2/96
C     Logic changed for EOM-CCSD gradients and added for EOM-MBPT(2) gradients
C     
      ELSEIF (IFLAGS2(h_IFLAGS2_program) .EQ. 3) THEN
         CCSD = .TRUE.
         EOM=.TRUE.
         write(6,3021)         
      ELSEIF (IFLAGS(h_IFLAGS_excite).GE.3) THEN
         CCSD=.TRUE.
         NAME3='CCSD'
         write(6,3020)NAME3
         EOM=.TRUE.
      ELSEIF(METHOD.EQ.1) THEN
         MBPT2=.TRUE.
         NAME='MBPT(2)'
         write(6,3000)NAME
      ELSE IF(METHOD.EQ.2) THEN
         MBPT3=.TRUE.
         NAME='MBPT(3)'
         write(6,3000)NAME
      ELSE IF(METHOD.EQ.3) THEN
         M4SDQ=.TRUE.
         NAME1='SDQ-MBPT(4)'
         write(6,3001) NAME1   
      ELSE IF(METHOD.EQ.4) THEN
         M4SDTQ=.TRUE.
         TRIP1=.TRUE.
         NAME='MBPT(4)'
         write(6,3000)NAME
      ELSE IF(METHOD.EQ.7) THEN
         UCC=.TRUE.
         NAME5='UCCSD(4)'
         write(6,3010)NAME5
      ELSE IF(METHOD.EQ.8 .OR. METHOD .EQ. 5) THEN
         CCD=.TRUE.
         NAME2='CCD'
         write(6,3002)NAME2
      ELSE IF(METHOD.EQ.9) THEN
         UCC=.TRUE.
         TRIP1=.TRUE.
         NAME7='UCC(4)'
         write(6,3012)NAME7
      ELSE IF(METHOD.EQ.10) THEN
         CCSD=.TRUE.
         NAME3='CCSD'
         write(6,3003)NAME3
      ELSE IF(METHOD.EQ.11) THEN
         CCSD=.TRUE.
         NAME6='CCSD+T(CCSDT)'
         TRIP1=.TRUE.
         write(6,3011) NAME6
      ELSE IF(METHOD.EQ.21) THEN
         QCISD=.TRUE.
         NAME5='QCISD(T)'
         TRIP1=.TRUE.
         TRIP2=.TRUE.
         write(6,3010) NAME5
      ELSE IF(METHOD.EQ.22) THEN
         CCSD=.TRUE.
         NAME='CCSD(T)'
         TRIP1=.TRUE.
         TRIP2=.TRUE.
         write(6,3000) NAME
      ELSE IF(METHOD.EQ.23 .OR. METHOD .EQ. 25) THEN
         QCISD=.TRUE.
         NAME4='QCISD'
         write(6,3004)NAME4
         Write(6,"(a,a)") "  WARNING! CISD or QCISD density matrices",
     +                    " are not fully tested!!!"
      ELSE IF(METHOD.EQ.42) THEN
         CCSD=.TRUE.
         NAME5='aCCSD(T)'
         TRIP1=.TRUE.
         TRIP2=.TRUE.
         LTRP=.TRUE.
         write(6,3010) NAME5
      ENDIF
C
C     
C     SET FLAGS FOR DERIVATIVE CALCULATION
C     
C     CALCULATE DENSITY MATRIX
C     
      DENS=.TRUE.
      RELAXED = (IFLAGS(h_IFLAGS_density).EQ.0)
C     
C     CALCULATE GRADIENTS WITh RESPECT TO NUCLEAR COORDINATES
C Note that derivatives with respect to nuclear perturbations require
C the I intermediate. For first-order property calculations, however,
C only the density matrices are required.
C Further clarification: Grad is true 1. Analytical gradient with geometry
C opt or vibrational frequency, exact vibrational frequencies with 
C analytic Hessians, single point analytical gradint calculations (I
C think this can now be done by only looking at grad_calc keyword only.
C     
CSSS      GRAD=(IFLAGS2(h_IFLAGS2_grad_calc).EQ.1  .AND.
CSSS     &      (IFLAGS2(h_IFLAGS2_geom_opt).ne.0  .OR.
CSSS     &       IFLAGS(h_IFLAGS_vib).EQ.3)        .OR. 
CSSS     &      (IFLAGS(h_IFLAGS_vib).EQ.1  .OR.
CSSS     &       IFLAGS(h_IFLAGS_deriv_lev) .GT.1 .OR.
CSSS     &       IFLAGS2(h_IFLAGS2_grad_calc).EQ.1))
C Commented and redone as follows (convoluted and do not work), 07/09
C Ajith Perera

      GRAD  = (IFLAGS2(h_IFLAGS2_grad_calc).EQ.1)  
      PROPS = (IFLAGS(h_IFLAGS_props)      .EQ.1) 
      NTORB = (IFLAGS2(h_IFLAGS2_pccd_reference).EQ.3)

      If (GRAD .AND. .NOT.RELAXED .AND. .NOT. PROPS) Then
C        This is trapped in joda and is only for redundancy.
         IF (.NOT. NTORB) THEN
         Write(6,"(2a)") "The unrelaxed density is not allowed",
     &                   "for this type of calculation."
         Call Errex
         ENDIF 
      Endif
      IF (.NOT. RELAXED) 
     &    Write(6,"(2a)") "  As requested response only density",
     &                    " matrices are computed."
C     
C Let's check whether the NONHF flag is set. If it is set, it
C could be the user or as a result of reference (ROHF or QRHF).
C If NONHF is set but the reference is RHF or UHF, then there are
C two possiblities. 1. user wants to use a non-diagonal Fock
C matrix (at least one of the vrt-vrt, occ-occ is not diagonal or
C occ-vrt blocks have non-zero elements or both). 2. user may
C set NONHF flag with or without knowing that he/she has a true
C non HF problem (see above for conditions for a non-HF problem).
C This is a problematic case, but the code should be able to
C generate a density or I intermediate that is identical to
C the standard RHF/UHF case.
C
      TRULY_NONHF = (IFLAGS(h_IFLAGS_ref).GE.2 .OR.
     &               IFLAGS(h_IFLAGS_qrhf_gen) .NE. 0)

C For QRHF gradients we need additional work. Also, let's set
C QRHF logical here. Note that this is a local flag (alias) in the 
C sense it does know anything about the other flags. Similarly
C set other the aliases, ROHF, SEMI (semicanonical), CANON 
C and (perturb orbitals). 
C     
#ifdef _DEBUG
      Print*, "The NONHF flags", iflags(h_IFLAGS_nonhf)
#endif 
      IF(IFLAGS(h_IFLAGS_qrhf_gen).NE.0)THEN
C     
C     JDW KKB stuff. Additions are in lower case.
C     
         QRHF=.TRUE.
         call getrec(20,'JOBARC','NDROPGEO',1,ndrgeo)
         if (ndrgeo.eq.0) then
            CALL QRHFSET
         else
            iflags(h_IFLAGS_nonhf)=0
         endif
C     JDW END
C     
         WRITE(6,3005)
C     
      ENDIF
      IF (.NOT. RELAXED) THEN
      WRITE(6,"(2a)") "  The QRHF set is not used for the response",
     &                " density matrices."
      ENDIF 
c OLD
c     IF(IFLAGS(h_IFLAGS_nonhf).NE.0.AND.IFLAGS(h_IFLAGS_ref).NE.2) THEN
c NEW
      IF(IFLAGS(h_IFLAGS_nonhf).EQ.1) THEN 
c I don't know why this was done this way: .AND.IFLAGS(h_IFLAGS_ref).NE.2
c END
C     
         NONHF=.TRUE.
C     
         WRITE(6,3006)
C     
      ELSE
         If (TRULY_NONHF) Then
            Write(6,*) "This is a non-HF calculation, but NONHF=FALSE!"
            Call Errex
         Endif
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
C     CHECK IF PERTURBED CANONICAL ORBITALS ARE REQUESTED
C     
      IF(IFLAGS(h_IFLAGS_pert_orb).EQ.1) THEN
         CANON=.TRUE.
         WRITE(6,3009)
      ENDIF
C     
C     CHECK WHETHER G(AB,CD) SHOULD BE CALCULATED OR NOT
C     
      IF(IFLAGS(h_IFLAGS_gamma_abcd).EQ.1) THEN
         GABCD=.TRUE.
         WRITE(6,3014)
      ENDIF

      RETURN
      END
