      SUBROUTINE RCC_DRIVER(WORK,MAXCOR,IUHF)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION WORK(MAXCOR)
      DIMENSION ECORR(3)
      INTEGER RELCYC
     
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,QCISD,
     +        UCC,CC2,RCCD,DRCCD
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

C Note at the first-iteration, T1 vector is zero and T2 is MBPT(2).

      CALL RCC_SETMET 
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

C For UHF we have to work with Coulomb only integrals and symmetric 
C amplitudes. The entire ACES II structure built on antisymetric 
C integals and amplitudes break down. So, AAAA and BBBB lists 
C and routines that supports operation on them have been readjusted 
C to work with symmetric integrals and amplitudes. This of course
C matters only for UHF. RHF calcs (as long as they can be solely 
C determined by ABAB combination) can proceed normally. 

      CALL RCC_SETLST(WORK(INEXT),MAXCOR,IUHF)

#ifdef _DEBUG_LVL0
      CALL RCC_DBG_SETLST(WORK(INEXT),MAXCOR,IUHF)
      CALL RCC_DBG_EXPANDED_ABIJ(WORK(INEXT),MAXCOR,IUHF)
#endif 
      IF (IUHF .NE. 0) THEN 
         CALL RCC_DRIVE_D2T2(WORK(INEXT),MAXCOR,IUHF)
         CALL RCC_RNABIJ(WORK(INEXT),MAXCOR,IUHF,'T')
      ELSE
         IF (RCCD) THEN
            NSIZE = IDSYMSZ(1,ISYTYP(1,46),ISYTYP(2,46))
            CALL GETALL(WORK(INEXT),NSIZE,1,46)
            CALL PUTALL(WORK(INEXT),NSIZE,1,44)
            CALL RCC_ABIJ_4RHF(WORK(INEXT),MAXCOR,IUHF)
         ENDIF 
         CALL RCC_RNABIJ(WORK(INEXT),MAXCOR,IUHF,'T')
      ENDIF 

#ifdef _DEBUG_LVL0
      CALL RCC_AMPSUM(WORK(1),MAXCOR,IUHF,0,.FALSE.,'T')
CSSS      call check_t2(work(1),maxcor,iuhf)
#endif 
CSSS      CALL RCC_SETLST(WORK(INEXT),MAXCOR,IUHF)

      CALL RCC_CMPENG(WORK(INEXT),MAXCOR,43,0,ECORR,ENERGY(1,1),
     &                ENERGY(1,2),IUHF,1)
      ELAST = ENERGY(1,2)

      IF (MOD(IFLAGS(h_IFLAGS_cc_extrapol),2).EQ.0) THEN
         CALL DRRLE(WORK(INEXT),MAXCOR,IUHF,RLECYC,.FALSE.)
      END IF
      IF (IFLAGS(h_IFLAGS_cc_extrapol) .EQ.1) THEN
         CALL DIISLST(1,IUHF,.FALSE.)
      END IF

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

C Construct W(mb,ej)= <mb||ej> -1/2T2(jn,fb)<mn||ef> intermediate

      CALL RCC_GENINT(WORK(INEXT),MAXCOR,IUHF)
C
C Construct +P(ij)P(ab) Sum_me T2(im,ae) W(mb,ej) (DRRNG)
C
      CALL RCC_DORNG(WORK(INEXT),MAXCOR,IUHF)
C
C Denominator weigh the T2 increments. 
C
#ifdef _DEBUG_LVL0
      CALL RCC_DBG(WORK(INEXT),MAXCOR,IUHF,ICYCLE)
#endif 
      CALL RCC_NEWT2(WORK(INEXT),MAXCOR,IUHF)

      IF (DAMP_PARAMETER .NE. 1.0D0) THEN
         CALL DAMP_CC_RESIDUAL(WORK(INEXT),MAXCOR,IUHF,.FALSE.,
     &                         44,61,90,90,DAMP_PARAMETER)
      ENDIF
C
#ifdef _DEBUG_LVL0
      CALL DRTSTS(WORK(INEXT),MAXCOR,ICYCLE,IUHF,ICONVG,ICONTL,
     &            .FALSE.,0,'T')
#endif
      CALL RCC_CMPENG(WORK(INEXT),MAXCOR,60,2,ECORR,
     &                ENERGY(ICYCLE+1,1),ENERGY(ICYCLE+1,2),IUHF,1)

      ELAST = ENERGY(ICYCLE+1,2)
C
      IF (IFLAGS(h_IFLAGS_cc_extrapol) .EQ.1)  THEN
         CALL RCC_DODIIS0(WORK(INEXT),MAXCOR,IUHF,1,ICYCLE,
     &                    ICONVG,ICONTL,.FALSE.,44,61,90,0,90,2,70,
     &                    '     ',DAMP_PARAMETER)
      END IF

      IF (ICONVG.NE.0) THEN

         CALL DRMOVE(WORK(INEXT),MAXCOR,IUHF,0,.FALSE.)
         CALL RCC_RNABIJ(WORK(INEXT),MAXCOR,IUHF,'T')

      END IF
#ifdef _DEBUG_LVLM
      if (iuhf .ne.0)Then
      call check_rhf_rcc(Work(Inext),Maxcor,Iuhf)
      endif 
#endif 

      IF (ICONVG.EQ.0) THEN

         IF (DRCCD) CALL RCC_ABTOAA(WORK(INEXT),MAXCOR,IUHF,61)
         CALL RCC_CMPENG(WORK(INEXT),MAXCOR,43,0,ECORR,
     &                   ENERGY(ICYCLE+1,1),ENERGY(ICYCLE+1,2),
     &                   IUHF,1)
     
         CALL RCC_AMPSUM(WORK(INEXT),MAXCOR,IUHF,0,.FALSE.,'T')
         IF (IUHF.NE.0) CALL S2PROJ(WORK(INEXT),MAXCOR,IUHF,.FALSE.)
C
         CALL RCC_FINISH(ICYCLE+1)
         CALL DDMPTGSS(WORK(INEXT), MAXCOR/IINTFP, IUHF, 0, 'TGUESS  ')
         CALL ACES_FIN 
         WRITE(6,1020)
 1020 FORMAT(/,77('-'),/,32X,'Exiting xvcc',/,77('-'),/)

         CALL ACES_EXIT(0)

      ENDIF 

      IF (MOD(IFLAGS(h_IFLAGS_cc_extrapol),2).EQ.0) THEN
         CALL DRRLE(WORK(INEXT),MAXCOR,IUHF,RLECYC,.FALSE.)
         CALL DRMOVE(WORK(INEXT),MAXCOR,IUHF,0,.FALSE.)
         CALL RCC_RNABIJ(WORK(INEXT),MAXCOR,IUHF,'T')
      END IF
C
      ICYCLE = ICYCLE +  1
      IF (ICYCLE .LE. NCYCLE) GOTO 100
       
      IF (ICONVG .NE. 0) THEN
         CALL RCC_AMPSUM(WORK(INEXT),MAXCOR,IUHF,0,.FALSE.,'T')
         Write(6,"(a,a)") "The RCCD or DRCCD equations did not converge",
     &                     " in alloted number of cyclces."
         call aces_exit(1)
      ENDIF 
C
      RETURN
      END
