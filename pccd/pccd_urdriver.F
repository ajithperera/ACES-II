#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */

      SUBROUTINE PCCD_URDRIVER(ICORE,MAXCOR,IUHF)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ICORE(MAXCOR)
      DIMENSION ECORR(3)
      INTEGER RELCYC

      LOGICAL NONHF
      LOGICAL PCCD,CCD,LCCD
      LOGICAL CIS,EOM
      LOGICAL BRUECK,BRKCNV
      LOGICAL DORESET,OO_CC
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,QCISD
      LOGICAL UCC,CC2
C
#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */
#include "machsp.com"
#include "syminf.com" /* for nirrep */
#include "sym.com"
#include "flags.com"
#include "flags2.com"
#include "bwcc.com"
#include "sympop.com"
#include "files.com"

      COMMON/CALC/PCCD,CCD,LCCD
      COMMON /ENERGY/ ENERGY(500,2)
      COMMON/EXCITE/EOM,CIS
      COMMON/NHFREF/NONHF
      COMMON/EXTRAPO/RLE
CSSS      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
CSSS     &                QCISD,UCC,CC2

      CALL ZERO(ENERGY,1000)
      CALL PCCD_SETMET()
      I0      = 1
      ICONTL  = IFLAGS(h_IFLAGS_cc_conv)
      NCYCLE  = IFLAGS(h_IFLAGS_cc_maxcyc)
      NKEEP   = IFLAGS(h_IFLAGS_cc_exporder)
      IEVERY  = IFLAGS(h_IFLAGS_tamp_sum)
      BRUECK  = (IFLAGS(h_IFLAGS_brueckner) .EQ. 1)
      METHOD  =  IFLAGS(h_IFLAGS_calclevel)
      DAMP_PARAMETER = IFLAGS2(h_IFLAGS2_damp_cc)*0.010D0
      SING1   = .FALSE.
      OO_CC   = (IFLAGS(h_IFLAGS_calc) .EQ. 52 .OR.
     &           IFLAGS(h_IFLAGS_calc) .EQ. 53 .OR.
     &           IFLAGS(h_IFLAGS_calc) .EQ. 54)
      LUOUT   = 6
      MOINTS  = 50
C
      CALL SETLST(ICORE(I0),MAXCOR,IUHF)
      CALL INSING(10,IUHF)
      ICYCLE = 1
      RLECYC = 0
      CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
      IF (IFLAGS(h_IFLAGS_cc_extrapol) .EQ.1) THEN
        CALL PCCD_DIISLST(1,IUHF,.FALSE.,"R")
      END IF
      CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')

      ICONVG = 1

      WRITE(6,801)
  801 FORMAT(T3,' Initial T amplitudes: ')
      CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
      WRITE(6,800)
  800 FORMAT(T3,' Correlation energies computed from initial ',
     &          'T amplitudes: ')

      CALL CMPENG(ICORE(I0),MAXCOR,43,0,ECORR,ENERGY(1,1),ENERGY(1,2),
     &            IUHF,1)
      ELAST = ENERGY(1,2)

 100  CONTINUE
      RLECYC = RLECYC + 1
C
C Start the T2 lists by adding <ab||ij>, This is actually a
C contribution to T2.
C
      CALL INITIN(ICORE(I0),MAXCOR,IUHF)

#ifdef _DEBUG_LVL0
       call check_t2(icore(i0),Maxcor,Iuhf) 
#endif 
      IF (METHOD.GE.6.AND.METHOD.NE.8) THEN
         CALL INITSN(ICORE(I0),MAXCOR,IUHF)
      END IF
c
C GENERATE W AND F INTERMEDIATES.
C
      IF (LINCC) THEN
         INTTYP = 3
      ELSE
         INTTYP = 2
      END IF
C 
C Generate W(mbej), W(mnij) W(mbej, F(ea), F(ij) and F(me)
C intermediates.  
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)") " The checksums of all the intms before genint"
      Write(6,"(a,2l)") " Flags singles and non-hf:",SING1,NONHF
      Call checkintms(icore(i0),maxcor,Iuhf,1)
      write(6,"(a)") " The T1/T2 residuals" 
      call check_t2(icore(i0),Maxcor,Iuhf) 
#endif 
C
      CALL GENINT(ICORE(I0),MAXCOR,IUHF,INTTYP,DO_HBAR_4LCCSD)
#ifdef _DEBUG_LVLM
      Write(6,"(a)") " The checksums of all the intms after genint"
      Write(6,*)
      Call checkintms(icore(i0),maxcor,Iuhf,1)
#endif 
C
C COMPUTE F INTERMEDIATE CONTRIBUTION TO T2 INCREMENT.
C
C Evaluate P(ab)Sum_e T(ij,ae){F(b,e)-1/2 sum_m T(m,b)F(m,e)}
C contribution to T2 and Sum_e T(i,e)F(a,e) contribution to T1.
C
      CALL FEACONT(ICORE(I0),MAXCOR,IUHF)
#ifdef _DEBUG_LVL0
       call check_t2(icore(i0),Maxcor,Iuhf)
#endif
C
C Evaluate -P(ij)Sum_e T(im,ab){F(mj)-1/2 sum_m T(j,e)F(m,e)}
C contribution to T2 and Sum_e T(m,a)F(m,i) contribution to T1.
C
      CALL FMICONT(ICORE(I0),MAXCOR,IUHF)
#ifdef _DEBUG_LVL0
       call check_t2(icore(i0),Maxcor,Iuhf)
#endif
C
      IF ((METHOD.GT.9.AND.SING1).OR.(METHOD.EQ.6.AND.SING1)) THEN

C Evaluate Sum_e T(im,ae)F(m,e) contribution to T1. 
C
       CALL FMECONT(ICORE(I0),MAXCOR,IUHF,1)
       IF (IUHF.NE.0) CALL FMECONT(ICORE(I0),MAXCOR,
     &                                  IUHF,2)
C
#ifdef _DEBUG_LVL0
       call check_t2(icore(i0),Maxcor,Iuhf)
#endif
C 
C Evaluate -P(ab)Sum_m T(m,a)<mb||ij> contribution to T2.
C 
         CALL T1INT2A(ICORE(I0),MAXCOR,IUHF)
C
C Evaluate +P(ij)Sum_m T(i,e)<ab||ej> contribution to T2.
C
         CALL T1INT2B(ICORE(I0),MAXCOR,IUHF)
C
C Evaluate -Sum_nf T(n,f)<na||if> contribution to T1.
C
         CALL T1INT1(ICORE(I0),MAXCOR,IUHF,1)
         IF (IUHF.NE.0) CALL T1INT1(ICORE(I0),MAXCOR,
     &                             IUHF,2)
C
#ifdef _DEBUG_LVL0
       call check_t2(icore(i0),Maxcor,Iuhf)
#endif
      END IF
C 
C DO W INTERMEDIATE CONTRIBUTION TO T2 EQUATION
C 
C Evaluate the following three contributions to T2:
C
C   -P(ij)P(ab) Sum_me T(i,e)T(m,a) <mb||ej> (T12INT2)
C   +P(ij)P(ab) Sum_me T(i,e)T(m,a) W(mb,ej) (DRRNG)
C   +1/2 Sum_mn Tau(mn,ab)W(mn,ij) + 
C    1/2 Sum_ef Tau(ij,ef)W(ab,ef)           (DRLAD)
C
#ifdef _DEBUG_LVL0
       call check_t2(icore(i0),Maxcor,Iuhf) 
#endif 
      CALL DRE3EN(ICORE(I0),MAXCOR,IUHF,0)
C 
#ifdef _DEBUG_LVLM
      write(6,"(a)") " The AAII residual block befor os"
      call check_res(icore(i0),maxcor,iuhf)
#endif 
C
C COMPUTE THE CONTRIBUTION OF DOUBLES TO T1 AND
C DENOMINATOR WEIGHT THE NEW INCREMENTS
C (ONLY CCSD AND QCISD METHODS)
C 
      IF ((METHOD.GE.9.OR.METHOD.EQ.6.OR.METHOD.EQ.7)
     &    .AND. .NOT. OO_CC) THEN
C
C Evaluate the following three contributions to T1:
C
C   -1/2 Sum_mef T(im,ef)<ma||ef>  (T2T1AA1, T2T1AB1)
C   -1/2 Sum_men T(nm,ei)<nm||ei>  (T2T1AA2, T2T1AB2)
C
#ifdef _DEBUG_LVL0
       call check_t2(icore(i0),Maxcor,Iuhf) 
#endif 
C Also, Do the T1 = T1/{f(i,i) - f(a,a)} to get a new T1
C
         CALL E4S(ICORE(I0),MAXCOR,IUHF,EDUMMY)
      END IF
#ifdef _DEBUG_LVL0
      CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
#endif 

C Also the call E4S for pCCD (OO-CCD) and Brueckner. This will do T2 into T1
C terms. Ajith Perera,10/2021. 

       IF (OO_CC .AND. BRUECK) THEN
          CALL E4S(ICORE(I0),MAXCOR,IUHF,EDUMMY)
       ENDIF 
C
C DENOMINATOR WEIGHT T2 INCREMENTS TO FORM NEW T2.
C 
C Do the T2 = T2/{f(i,i) + f(j,j) - f(a,a) - f(b,b)}
C to get a new T2.
C
      IF (pCCD) THEN
C Ajith Perera, 10/2021.
         Write(6,"(20x,a)") " ----------Warning------------"
         Write(6,"(2a)")  " A pCCD or pCCD like calculation is being",
     +                    " performed and the off-diagonal"
         Write(6,"(a)")   " blocks of T2 is set to zero"
         Write(6,"(20x,a)") " -----------------------------"
         call pccd_reset_vcc(icore(i0),Maxcor,Iuhf,63)
      ENDIF 
#ifdef _DEBUG_LVL0
       call check_t2(icore(i0),Maxcor,Iuhf)
#endif
      CALL NEWT2(ICORE(I0),MAXCOR,IUHF)

      CALL CMPENG(ICORE(I0),MAXCOR,60,2,ECORR,ENERGY(ICYCLE+1,1),
     &            ENERGY(ICYCLE+1,2),IUHF,1)
      ELAST = ENERGY(ICYCLE+1,2)
C
C IF ICONVG IS EQUAL ZERO, CONVERGENCE HAS BEEN ACHIEVED, EXIT
C VIA FINISH
C
      IF (IFLAGS(h_IFLAGS_cc_extrapol) .EQ.1) THEN
         CALL PCCD_DODIIS0(ICORE(I0),MAXCOR/IINTFP,IUHF,1,ICYCLE,
     &                ICONVG,ICONTL,.FALSE.,44,61,90,0,90,2,70)
      END IF

      IF (ICONVG.NE.0) THEN
         CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,0,SING1)
         CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
      END IF

      IF (ICONVG.EQ.0) THEN
         CALL CMPENG(ICORE(I0),MAXCOR,43,0,ECORR,ENERGY(ICYCLE+1,1),
     &               ENERGY(ICYCLE+1,2),IUHF,1)
         CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
C
CSSS         CALL FINISH(ICYCLE+1,PCCD,.FALSE.,.FALSE.)
         CALL PCCD_UFINISH(ICYCLE+1,"R")
         CALL INIT2(IUHF,PCCD.OR.CCD.OR.LCCD)
         CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,100,SING1)
         CALL DDMPTGSS(ICORE(I0), MAXCOR/IINTFP, IUHF, 0, 'TGUESS  ')
C
         IF (DORESET) THEN
            IF(IUHF.EQ.0)CALL RESET(ICORE(I0),MAXCOR,IUHF)
         END IF 
        WRITE(6,1020)
 1020 FORMAT(/,77('-'),/,32X,'Exiting pccd_urdriver',/,77('-'),/)

         RETURN
      ENDIF 
C
      ICYCLE = ICYCLE +  1
      IF (ICYCLE .LE. NCYCLE) GOTO 100

      IF (ICONVG .NE. 0) THEN
         CALL PCCD_AMPSUM(ICORE(I0),MAXCOR,IUHF,0,.FALSE.,'T')
         IF (pCCD) THEN
         Write(6,"(2a)") "The pCCD equations did not converge",
     &                     " in alloted number of cyclces."
         ELSE IF(CCD) THEN
         Write(6,"(2a)") "The CCD equations did not converge",
     &                     " in alloted number of cyclces."
         ELSE IF(LCCD) THEN
         Write(6,"(2a)") "The LCCD equations did not converge",
     &                     " in alloted number of cyclces."
         ENDIF
      ENDIF



      RETURN
      END

