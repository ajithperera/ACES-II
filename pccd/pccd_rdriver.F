      SUBROUTINE PCCD_RDRIVER(WORK,MAXCOR,IUHF)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION WORK(MAXCOR)
      DIMENSION ECORR(3)
      INTEGER RELCYC
     
      LOGICAL PCCD,CCD,LCCD
      LOGICAL NONHF
      LOGICAL CIS,EOM
      LOGICAL BRUECK,BRKCNV
C
#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */
#include "machsp.com"
#include "syminf.com" /* for nirrep */
#include "sym.com"
#include "flags.com"
#include "flags2.com"
#include "bwcc.com"
#include "sympop.com"

      COMMON/CALC/PCCD,CCD,LCCD
      COMMON /ENERGY/ ENERGY(500,2)
      COMMON/EXCITE/EOM,CIS
      COMMON/NHFREF/NONHF
      COMMON/EXTRAPO/RLE

C Note at the first-iteration, T1 vector is zero and T2 is MBPT(2).

      CALL PCCD_SETMET()
C
C Initializations of various intermediate lists etc.
C
      INEXT   = 1
      RLECYC  = 0
      IREF    = 1
      IBWPASS = 0
      BRUECK  = .FALSE.
      ICONTL  = IFLAGS(h_IFLAGS_cc_conv)
      NCYCLE  = IFLAGS(h_IFLAGS_cc_maxcyc)
      NKEEP   = IFLAGS(h_IFLAGS_cc_exporder)
      IEVERY  = IFLAGS(h_IFLAGS_tamp_sum)
      BRUECK  = (IFLAGS(h_IFLAGS_brueckner) .EQ. 1)

      CALL PCCD_SETLST(WORK(INEXT),MAXCOR,IUHF,"R")

      CALL PCCD_INIT(WORK(INEXT),MAXCOR,IUHF,"R")
      CALL PCCD_RNABIJ(WORK(INEXT),MAXCOR,IUHF,'T')

      IF (IFLAGS(h_IFLAGS_cc_extrapol) .EQ.1) THEN
         CALL PCCD_DIISLST(1,IUHF,.FALSE.,"R")
      END IF

      ICYCLE = 1

 100  CONTINUE
     
      RELCYC = RLECYC + 1
      CALL PCCD_INITIN(WORK(INEXT),MAXCOR,IUHF)
      CALL PCCD_INITSN(WORK(INEXT),MAXCOR,IUHF,ICYCLE,"T")

#ifdef _DEBUG_LVL0
       call pccd_check_t2(work(inext),Maxcor,Iuhf)
#endif

      CALL PCCD_CMPENG(WORK(INEXT),MAXCOR,46,ECORR,ENERGY(1,1))

C Construct W(mb,ej)= <mb||ej> -1/2T2(jn,fb)<mn||ef> intermediate

      CALL PCCD_GENINT(WORK(INEXT),MAXCOR,IUHF)
C
C Evaluate P(ab)Sum_e T(ij,ae){F(b,e)-1/2 sum_m T(m,b)F(m,e)}
C contribution to T2 and Sum_e T(i,e)F(a,e) contribution to T1.
C
      CALL PCCD_FEACONT(WORK(INEXT),MAXCOR,IUHF)

#ifdef _DEBUG_LVL0
       call pccd_check_t2(work(inext),Maxcor,Iuhf)
#endif
C
C Evaluate -P(ij)Sum_e T(im,ab){F(mj)-1/2 sum_m T(j,e)F(m,e)}
C contribution to T2 and Sum_e T(m,a)F(m,i) contribution to T1.
C
      CALL PCCD_FMICONT(WORK(INEXT),MAXCOR,IUHF)

#ifdef _DEBUG_LVL0
       call pccd_check_t2(work(inext),Maxcor,Iuhf)
#endif
C
C Evaluate Sum_e T(im,ae)F(m,e) contribution to T1.
C  
#ifdef _T1ON
      IF (IFLAGS(h_IFLAGS_calc) .GT. 9 .AND.S) THEN
          CALL PCCD_FMECONT(WORK(INEXT),MAXCOR,IUHF,1)
          IF (IUHF.NE.0) CALL PCCD_FMECONT(WORK(INEXT),
     &                                 MAXCOR,IUHF,2)
      ENDIF 
#endif


C   +P(ij)P(ab) Sum_me T(i,e)T(m,a) W(mb,ej) (DRRNG)
C   +1/2 Sum_mn Tau(mn,ab)W(mn,ij) +
C    1/2 Sum_ef Tau(ij,ef)W(ab,ef)           (DRLAD)
C
      CALL PCCD_DRE3EN(WORK(INEXT),MAXCOR,IUHF,0)

C Evaluate the following three contributions to T1:
C
C   -1/2 Sum_mef T(im,ef)<ma||ef>  (T2T1AA1, T2T1AB1)
C   -1/2 Sum_men T(nm,ei)<nm||ei>  (T2T1AA2, T2T1AB2)
C
#ifdef _DEBUG_LVL0
       call pccd_check_t2(work(inext),Maxcor,Iuhf)
#endif
C Also, Do the T1 = T1/{f(i,i) - f(a,a)} to get a new T1
#ifdef _T1ON
       IF (NONHF) CALL PCCD_E4S(WORK(INEXT),MAXCOR,IUHF,
     &                          EDUMMY)
#endif 
C
C Denominator weigh the T2 increments. 
C
#ifdef _HF_ONLY
      CALL PCCD_SET2NULL(WORK(INEXT),MAXCOR,IUHF,63)
#endif 
      IF(PCCD) CALL PCCD_RESET(WORK(INEXT),MAXCOR,IUHF,63)
      CALL PCCD_NEWT2(WORK(INEXT),MAXCOR,IUHF)

#ifdef _DEBUG_LVL0
       call pccd_check_t2(work(inext),Maxcor,Iuhf)
#endif
C
#ifdef _DEBUG_LVL0
      CALL PCCD_DRTSTS(WORK(INEXT),MAXCOR,ICYCLE,IUHF,ICONVG,ICONTL,
     &            .FALSE.,0,'T')
#endif

      IF (IFLAGS(h_IFLAGS_cc_extrapol) .EQ.1)  THEN
         CALL PCCD_DODIIS0(WORK(INEXT),MAXCOR,IUHF,1,ICYCLE,
     &                    ICONVG,ICONTL,.FALSE.,44,61,90,0,90,2,70)

      END IF
      CALL PCCD_CMPENG(WORK(INEXT),MAXCOR,63,ECORR,ENERGY(ICYCLE+1,1))

      IF (ICONVG.NE.0) THEN

         CALL PCCD_DRMOVE(WORK(INEXT),MAXCOR,IUHF,0)
         CALL PCCD_RNABIJ(WORK(INEXT),MAXCOR,IUHF,'T')

      END IF

      IF (BRUECK .AND. (ICYCLE.GT.1)) THEN
         CALL CHKBKSTP(WORK(INEXT),MAXCOR/IINTFP,IUHF,ICONVG,pCCD)
      END IF

      IF (ICONVG.EQ.0) THEN
         CALL PCCD_INITL2D2(IUHF)
         CALL PCCD_DRMOVE(WORK(INEXT),MAXCOR,IUHF,100)
         CALL PCCD_CMPENG(WORK(INEXT),MAXCOR,46,ECORR,ENERGY(1,1))
         CALL PCCD_AMPSUM(WORK(INEXT),MAXCOR,IUHF,0,.FALSE.,'T')

         CALL PCCD_RFINISH(ICYCLE+1,"R")
         CALL PCCD_DDMPTGSS(WORK(INEXT), MAXCOR/IINTFP, IUHF, 0,
     &                      'TGUESS  ')
         IF (BRUECK) CALL BRUECKIT(WORK(INEXT),MAXCOR,IUHF,
     &                             BRKCNV,pCCD)
         WRITE(6,1020)
 1020 FORMAT(/,77('-'),/,32X,'Exiting pccd_rdriver',/,77('-'),/)
         RETURN
      ENDIF 
C
      ICYCLE = ICYCLE +  1
      IF (ICYCLE .LE. NCYCLE) GOTO 100
       
      IF (ICONVG .NE. 0) THEN
         CALL PCCD_AMPSUM(WORK(INEXT),MAXCOR,IUHF,0,.FALSE.,'T')
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
C
      RETURN
      END
