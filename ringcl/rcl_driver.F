      SUBROUTINE RCL_DRIVER(WORK,MAXCOR,IUHF,SOLVE_4LAMBDA)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION WORK(MAXCOR)
      DIMENSION ECORR(3)
     
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,QCISD,
     +        UCC,CC2,RCCD,DRCCD
      LOGICAL SOLVE_4LAMBDA
C
C A ring CCD (rCCD) and direct ring CCD (drCCD) as described in JCP 139, 
C 104113 (2013) is implemented. 
C
#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */
#include "machsp.com"
#include "syminf.com" /* for nirrep */
#include "sym.com"
#include "flags.com"
#include "flags2.com"
#include "bwcc.com"
#include "sympop.com"

      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
      COMMON /ENERGY/ ENERGY(500,2),IXTRLE(500)

C If only RCC and DRCC Hbar is required then do not need to 
C solve for lambda equations. There are instances we need both (example
C RCC and DRCC property calculations.

      IF ( .NOT. SOLVE_4LAMBDA) RETURN

C Note at the first-iteration, T1 vector is zero and T2 is MBPT(2).

      CALL RCL_SETMET(.FALSE.)
C
C Initializations of various intermediate lists etc.
C
      INEXT   = 1
      RLECYC  = 0
      IREF    = 1
      IBWPASS = 0
      ICONTL  = IFLAGS(h_IFLAGS_cc_conv)
      NCYCLE  = IFLAGS(h_IFLAGS_cc_maxcyc)
      NKEEP   = IFLAGS(h_IFLAGS_cc_exporder)
      IEVERY  = IFLAGS(h_IFLAGS_tamp_sum)
      IRREPX  = 1
C For UHF we have to work with Coulomb only integrals and symmetric 
C amplitudes. The entire ACES II structure built on antisymetric 
C integals and amplitudes break down. So, AAAA and BBBB lists 
C and routines that supports operation on them have been readjusted 
C to work with symmetric integrals and amplitudes. This of course
C matters only for UHF. RHF calcs (as long as they can be solely 
C determined by ABAB combination) can proceed normally. 

      CALL RCL_SETLST(WORK(INEXT),MAXCOR,IUHF)
      
      NSIZE = IDSYMSZ(1,ISYTYP(1,46),ISYTYP(2,46))
      CALL GETALL(WORK(INEXT),NSIZE,IRREPX,61)
      CALL PUTALL(WORK(INEXT),NSIZE,IRREPX,144)

      IF (IUHF .NE. 0) THEN
         CALL GETALL(WORK(INEXT),NSIZE,IRREPX,62)
         CALL PUTALL(WORK(INEXT),NSIZE,IRREPX,145)
      ENDIF

      CALL GETALL(WORK(INEXT),NSIZE,IRREPX,63)
      CALL PUTALL(WORK(INEXT),NSIZE,IRREPX,146)

      CALL RCC_RNABIJ(WORK(INEXT),MAXCOR,IUHF,'L')

#ifdef _DEBUG_LVL0
      CALL RCL_CHECK_UHF(WORK(1),MAXCOR,IUHF)
      CALL RCC_AMPSUM(WORK(1),MAXCOR,IUHF,100,.FALSE.,'L')
#endif 
C Strictly speaking the following contribtion must not be there.
C See the calls rcl_fabinl2 and rcl_fijinl2 below how they 
C contribute to lambda. The W(a<b;i<j) expanded is constructed
C and saved before we compute the Hbar(a,b) and Hbar(i,j) 
C contributions.

      SCALE = 1.0D0
      Call RCL_EXPANDED_ABIJ(WORK(INEXT),MAXCOR,IUHF,SCALE)

C Construct Hbar(a,b) = f(a,b) - 1/2T2(mn,ae)<mn||be> (f(a,b) is
C part of the denominator arrays and not included).

      CALL RCL_FORM_HBAR_AB(WORK(INEXT),MAXCOR,IUHF,SCALE)

C Construct Hbar(j,i) = f(j,i) + 1/2T2(im,ef)<jm||ef> (f(i,j) is
C part of the denominator arrays and not included).

      
      CALL RCL_FORM_HBAR_JI(WORK(INEXT),MAXCOR,IUHF,Scale)

      CALL RCC_CMPENG(WORK(INEXT),MAXCOR,143,0,ECORR,ENERGY(1,1),
     &                ENERGY(1,2),IUHF,1)
      CALL RCC_CMPENG(WORK(INEXT),MAXCOR,60,0,ECORR,ENERGY(1,1),
     &                ENERGY(1,2),IUHF,1)

      ELAST = ENERGY(1,2)

      IF (MOD(IFLAGS(h_IFLAGS_cc_extrapol),2).EQ.0) THEN
         CALL DRRLE(WORK(INEXT),MAXCOR,IUHF,RLECYC,.TRUE.)
      END IF
      IF (IFLAGS(h_IFLAGS_cc_extrapol) .EQ.1) THEN
         CALL DIISLST(1,IUHF,.FALSE.)
      END IF
C           - 
C Construct W(ej,mb)= W(ej,mb) + 1/2T2(mn,ef)<nj||bf> intermediate

      IF (DRCCD) THEN
         CALL DRCL_DWMBEJ(WORK(INEXT),MAXCOR,IUHF)
      ELSEIF (RCCD) THEN
         IF (IUHF.EQ.0) THEN
            CALL RCL_DWMBEJ_R(WORK(INEXT),MAXCOR,IUHF)
         ELSE 
            CALL RCL_DWMBEJ_U(WORK(INEXT),MAXCOR,IUHF)
         ENDIF 
      ENDIF 

      ICYCLE = 1
      DAMP_PARAMETER = IFLAGS2(h_IFLAGS2_damp_cc)*0.010D0

#ifdef _DEBUG_LVL0
      Write(*,*)
      write(*,"(3x,a,F6.4)") "The RCC damp parameter = ", DAMP_PARAMETER
      Write(*,*)
#endif

 100  CONTINUE
     
      RELCYC = RLECYC + 1

      IF (DRCCD) THEN
         CALL DRCC_INITIN(WORK(INEXT),MAXCOR,IUHF)
      ELSE IF (RCCD) THEN
         CALL RCC_INITIN(WORK(INEXT),MAXCOR,IUHF)
      ENDIF 

#ifdef _DEBUG_LVLM
      call check_l2(work(inext),maxcor,iuhf)
#endif 
C These contributions are not in the strict formulation of lambda, but
C would be if we assume a e^T form for ring-CCD wave function.

#ifdef _EXPONENTIAL_FORM
     CALL RCL_FAB_INL2(WORK(INEXT),MAXCOR,IUHF)
     CALL RCL_FJI_INL2(WORK(INEXT),MAXCOR,IUHF)
#endif 

#ifdef _DEBUG_LVLM
      CALL RCC_CMPENG(WORK(INEXT),MAXCOR,60,0,ECORR,ENERGY(1,1),
     &                ENERGY(1,2),IUHF,1)
      call check_l2(work(inext),maxcor,iuhf)
      if (iuhf.eq.0)Then
      call check_rhf_rcc(Work(Inext),Maxcor,Iuhf)
      endif
#endif 
C                                        -
C Construct L2(ij,ab)=P(ij)P(ab)L2(ae,im)W(ej,mb)

      CALL RCL_DORNG(WORK(INEXT),MAXCOR,IUHF)
C
C Denominator weigh the L2 increments. 
C
      CALL RCL_NEWT2(WORK(INEXT),MAXCOR,IUHF)

      IF (DAMP_PARAMETER .NE. 1.0D0) THEN
         CALL DAMP_CC_RESIDUAL(WORK(INEXT),MAXCOR,IUHF,.FALSE.,
     &                         144,61,190,190,DAMP_PARAMETER)
      ENDIF
C
C Compute RCC correlation energies for each cycle.
C
#ifdef _DEBUG_LVL0
      CALL RCL_CHECK_UHF(WORK(1),MAXCOR,IUHF)
      CALL DRTSTS(WORK(INEXT),MAXCOR,ICYCLE,IUHF,ICONVG,ICONTL,
     &            .FALSE.,100,'L')
#endif
      CALL RCC_CMPENG(WORK(INEXT),MAXCOR,60,2,ECORR,
     &                ENERGY(ICYCLE+1,1),ENERGY(ICYCLE+1,2),IUHF,1)

      ELAST = ENERGY(ICYCLE+1,2)
C
      IF (IFLAGS(h_IFLAGS_cc_extrapol) .EQ.1)  THEN
         CALL RCC_DODIIS0(WORK(INEXT),MAXCOR,IUHF,1,ICYCLE,
     &                    ICONVG,ICONTL,.FALSE.,144,61,190,0,190,2,70,
     &                    '     ',DAMP_PARAMETER)
      ENDIF 

      IF (ICONVG .NE. 0) THEN
         CALL DRMOVE(WORK(INEXT),MAXCOR,IUHF,100,.FALSE.)
         CALL RCC_RNABIJ(WORK(INEXT),MAXCOR,IUHF,'L')
      END IF

#ifdef _DEBUG_LVLM
      if (iuhf .ne.0)Then
      call check_rhf_rcc(Work(Inext),Maxcor,Iuhf)
      endif
#endif

      IF (ICONVG.EQ.0) THEN
         IF (DRCCD) CALL RCC_ABTOAA(WORK(INEXT),MAXCOR,IUHF,61)
         CALL RCC_CMPENG(WORK(INEXT),MAXCOR,143,0,ECORR,
     &                   ENERGY(ICYCLE+1,1),ENERGY(ICYCLE+1,2),
     &                   IUHF,1)
     
         CALL RCC_AMPSUM(WORK(INEXT),MAXCOR,IUHF,0,.FALSE.,'L')
C
         CALL RCC_FINISH(ICYCLE+1)
         CALL ACES_FIN 
         WRITE(6,1020)
 1020 FORMAT(/,77('-'),/,32X,'Exiting xvcc',/,77('-'),/)

         CALL ACES_EXIT(0)

      ENDIF 

      IF (MOD(IFLAGS(h_IFLAGS_cc_extrapol),2).EQ.0) THEN
         CALL DRRLE(WORK(INEXT),MAXCOR,IUHF,RLECYC,.FALSE.)
         CALL DRMOVE(WORK(INEXT),MAXCOR,IUHF,100,.FALSE.)
         CALL RCC_RNABIJ(WORK(INEXT),MAXCOR,IUHF,'L')
      END IF
C
      ICYCLE = ICYCLE +  1
      IF (ICYCLE .LE. NCYCLE) GOTO 100
       
      IF (ICONVG .NE. 0) THEN
         CALL RCC_AMPSUM(WORK(INEXT),MAXCOR,IUHF,0,.FALSE.,'L')
         Write(6,"(2a)") "The RCCD or DRCCD equations did not converge",
     &                   " in alloted number of cyclces."
         call aces_exit(1)
      ENDIF 
       
C
      RETURN
      END
