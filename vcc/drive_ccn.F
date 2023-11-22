C
      SUBROUTINE DRIVE_CCN(ICORE,MAXCOR,IUHF,NCYCLE,ICONTL,
     &                     DAMP_PARAMETER)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION ICORE(MAXCOR)
      DIMENSION ECORR(3)
     
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC,CC2
      LOGICAL DO_HBAR_4LCCSD
C
C A CC2 described as in CPL,243,409-418 (1995) is implemented. 
C
#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */
#include "machsp.com"
#include "syminf.com" /* for nirrep */
#include "sym.com"
#include "jodaflags.com"

      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC,CC2
      COMMON /ENERGY/ ENERGY(500,2),IXTRLE(500)
C
C Note at the first-iteration, T1 vector is zero and T2 is MBPT(2).
C
#ifdef _DEBUG_LVL0
      CALL AMPSUM(ICORE(1),MAXCOR,IUHF,0,.TRUE.,'T')
#endif 
C
C Initializations of various intermediate lists etc.
C
      INEXT = 1
      SING1 = .FALSE.
      DO_HBAR_4LCCSD = .FALSE.

      CALL RNABIJ(ICORE(INEXT),MAXCOR,IUHF,'T')
      CALL SETLST(ICORE(INEXT),MAXCOR,IUHF)

      IF (MOD(IFLAGS(21),2).EQ.0) THEN
         CALL DRRLE(ICORE(I0),MAXCOR,IUHF,RLECYC,.FALSE.)
      END IF
      IF (IFLAGS(21).EQ.1) THEN
         CALL DIISLST(1,IUHF,.TRUE.)
      END IF

      RELCYC = 0
      ICYCLE = 1

 100  CONTINUE
     
      RELCYC = RLECYC + 1

      CALL INITIN(ICORE(INEXT),MAXCOR,IUHF)
      CALL INITSN(ICORE(INEXT),MAXCOR,IUHF)
C
C Generate W(mbej), W(mnij) W(mbej, F(ea), F(ij) and F(me)
C intermediates. CC2 like methods requires that the next to last 
C argument be 2. 
C
      CALL GENINT(ICORE(INEXT),MAXCOR,IUHF,2,DO_HBAR_4LCCSD)
C 
C Evalaute Sum_e T(i,e)F(a,e) contribution to T1.
C
      CALL FEACONT_MODF(ICORE(INEXT),MAXCOR,IUHF)
C
C Evaluate Sum_e T(m,a)F(m,i) contribution to T1.
C
      CALL FMICONT_MODF(ICORE(INEXT),MAXCOR,IUHF)
C
C Evalaute Sum_e T(m,a)F(m,i) contribution to T1.
C  
      IF (SING1) THEN

         CALL FMECONT(ICORE(INEXT),MAXCOR,IUHF,1)
         IF(IUHF.NE.0) CALL FMECONT(ICORE(INEXT),MAXCOR,IUHF,2)
C
C Evaluate -P(ab)Sum_m T(m,a)<mb||ij> contribution to T2.
C
         CALL T1INT2A(ICORE(INEXT),MAXCOR,IUHF)
C
C Evaluate +P(ij)Sum_m T(i,e)<ab||ej> contribution to T2.
C
         CALL T1INT2B(ICORE(INEXT),MAXCOR,IUHF)
C
C Evaluate -Sum_nf T(n,f)<na||if> contribution to T1.
C
         CALL T1INT1(ICORE(INEXT),MAXCOR,IUHF,1)
         IF(IUHF.NE.0) CALL T1INT1(ICORE(INEXT),MAXCOR,IUHF,2)

      ENDIF
C
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
      CALL DRE3EN(ICORE(INEXT),MAXCOR,IUHF,0)
C
C Evaluate the following three contributions to T1: Also,
C denominator weigh the T1 increments
C
C   -1/2 Sum_mef T(im,ef)<ma||ef>  (T2T1AA1, T2T1AB1)
C   -1/2 Sum_men T(nm,ei)<nm||ei>  (T2T1AA2, T2T1AB2)
C
      CALL E4S(ICORE(INEXT),MAXCOR,IUHF,DUMMY)
C 
C Denominator weigh the T2 increments. 
C
      CALL NEWT2(ICORE(INEXT),MAXCOR,IUHF)
C
C At this point, we should have all the t1 contributions 
C computed. So, we can look at the difference between
C the current and previouse T1 vector. 
C 
#ifdef _DEBUG_LVL0
      CALL AMPSUM(ICORE(INEXT),MAXCOR,IUHF,0,.TRUE.,'T')
      CALL DRTSTS(ICORE(INEXT),MAXCOR,ICYCLE,IUHF,ICONVG,ICONTL,
     &            SING1,0,'T')
#endif 
C
C Compute correlation energies for each cycle.
C
      CALL CMPENG(ICORE(INEXT),MAXCOR,60,2,ECORR,ENERGY(ICYCLE+1,1),
     &            ENERGY(ICYCLE+1,2),IUHF,1)
C

      IF (IFLAGS(21).EQ.1) THEN
         CALL DODIIS0(ICORE(INEXT),MAXCOR/IINTFP,IUHF,1,ICYCLE,
     &                ICONVG,ICONTL,SING1,44,61,90,0,90,2,70,
     &                '     ',DAMP_PARAMETER)
      END IF

      IF (ICONVG.NE.0) THEN

         CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,0,SING1)
         CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')

      END IF

      IF (ICONVG.EQ.0) THEN

         CALL CMPENG(ICORE(INEXT),MAXCOR,43,0,ECORR,ENERGY(ICYCLE+1,1),
     &               ENERGY(ICYCLE+1,2),IUHF,1)
     
         CALL AMPSUM(ICORE(INEXT),MAXCOR,IUHF,0,SING1,'T')
         IF (IUHF.NE.0) CALL S2PROJ(ICORE(INEXT),MAXCOR,IUHF,SING1)
C
         CALL FINISH(ICYCLE+1)
         CALL DDMPTGSS(ICORE(INEXT), MAXCOR/IINTFP, IUHF, 0, 'TGUESS  ')

         CALL INIT2(IUHF)
         CALL DRMOVE(ICORE(INEXT),MAXCOR,IUHF,100,SING1)
         IF (IUHF.EQ.0.AND..NOT.UCC) CALL RESET(ICORE(INEXT),
     &                                          MAXCOR,IUHF)

         CALL ACES_FIN 
         WRITE(6,1020)
 1020 FORMAT(/,77('-'),/,32X,'Exiting xvcc',/,77('-'),/)

         CALL ACES_EXIT(0)

      ENDIF 

      IF (MOD(IFLAGS(21),2).EQ.0) THEN
         CALL DRRLE(ICORE(INEXT),MAXCOR,IUHF,RLECYC,.FALSE.)
         CALL DRMOVE(ICORE(INEXT),MAXCOR,IUHF,0,SING1)
         CALL RNABIJ(ICORE(INEXT),MAXCOR,IUHF,'T')
      END IF
C
      ICYCLE = ICYCLE +  1
      IF (ICYCLE .LE. NCYCLE) GOTO 100
       
      IF (ICONVG .NE. 0) THEN
         CALL AMPSUM(ICORE(INEXT),MAXCOR,IUHF,0,SING1,'T')
         Write(6,"(a,a)") "The CC2 equations did not converge",
     &                     " in alloted number of cyclces."
         call aces_exit(1)
      ELSE
 
         Write(6,"(A,1x,i3,1x,A)")"T1 iterations converged in",
     &                             ICYCLE, "cyclces."
      ENDIF 
       
C
      RETURN
      END
