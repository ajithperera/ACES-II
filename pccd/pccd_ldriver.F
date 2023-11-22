      SUBROUTINE PCCD_LDRIVER(WORK,MAXCOR,IUHF)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION WORK(MAXCOR)
      DIMENSION ECORR(3)
      INTEGER RELCYC
     
      LOGICAL PCCD,CCD,LCCD
      LOGICAL NONHF
      LOGICAL CIS,EOM
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

C Initializations of various intermediate lists etc.
C
      INEXT   = 1
      RLECYC  = 0
      IREF    = 1
      IBWPASS = 0
      IRREPX  = 1
      ICONTL  = IFLAGS(h_IFLAGS_cc_conv)
      NCYCLE  = IFLAGS(h_IFLAGS_cc_maxcyc)
      NKEEP   = IFLAGS(h_IFLAGS_cc_exporder)
      IEVERY  = IFLAGS(h_IFLAGS_tamp_sum)

      CALL PCCD_SETLST(WORK(INEXT),MAXCOR,IUHF,"L")
      CALL PCCD_INIT(WORK(INEXT),MAXCOR,IUHF,"L")
      CALL PCCD_RNABIJ(WORK(INEXT),MAXCOR,IUHF,'L')

      IF (LCCD) THEN
         Write(6,"(a)") "      Lambda is constructed for LCCD."
         CALL MVSYM(WORK(INEXT),MAXCOR,IUHF,IRREPX,46,146)
         CALL GETLST(WORK(INEXT),1,1,1,1,90)
         CALL PUTLST(WORK(INEXT),1,1,1,1,190)
         CALL LGENINT(WORK(INEXT),MAXCOR,IUHF)
         WRITE(6,1020)
         RETURN
      ENDIF
     
      IF (IFLAGS(h_IFLAGS_cc_extrapol) .EQ.1) THEN
         CALL PCCD_DIISLST(1,IUHF,.FALSE.,"L")
      END IF
      CALL PCCD_FORMWL(WORK(INEXT),MAXCOR,IUHF)
      CALL PCCD_WTWTW(WORK(INEXT),MAXCOR,IUHF)
      CALL PCCD_FORMW4(WORK(INEXT),MAXCOR,IUHF,.TRUE.,.TRUE.,
     &                 .TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.,
     &                 .TRUE.)
      CALL PCCD_FORMW5(WORK(INEXT),MAXCOR,IUHF,.TRUE.,.TRUE.,
     &                 .TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.,
     &                 .TRUE.)
      ICYCLE = 1
 100  CONTINUE

      RELCYC = RLECYC + 1
      CALL PCCD_INITIN(WORK(INEXT),MAXCOR,IUHF)
      CALL PCCD_INITSN(WORK(INEXT),MAXCOR,IUHF,ICYCLE,"L")

      CALL PCCD_CMPENG(WORK(INEXT),MAXCOR,146,ECORR,ENERGY(1,1))

      CALL PCCD_LGENINT(WORK(INEXT),MAXCOR,IUHF)
C
#ifdef _DEBUG_LVL0
       call pccd_checkintms(work(inext),maxcor,iuhf,0)
       call pccd_check_t2(work(inext),Maxcor,Iuhf)
#endif
C
C Evaluate P(ab)Sum_e T(ij,ae){F(b,e)-1/2 sum_m T(m,b)F(m,e)}
C contribution to T2 and Sum_e T(i,e)F(a,e) contribution to T1.
C
      CALL PCCD_F1INL2(WORK(INEXT),MAXCOR,IUHF)
      CALL PCCD_F2INL2(WORK(INEXT),MAXCOR,IUHF)

#ifdef _DEBUG_LVL0
       call pccd_check_t2(work(inext),Maxcor,Iuhf)
#endif

      CALL PCCD_G1INL2(WORK(INEXT),MAXCOR,IUHF)
      CALL PCCD_G2INL2(WORK(INEXT),MAXCOR,IUHF)

#ifdef _DEBUG_LVL0
       call pccd_check_t2(work(inext),Maxcor,Iuhf)
#endif

C   +P(ij)P(ab) Sum_me T(i,e)T(m,a) W(mb,ej) (DRRNG)
C   +1/2 Sum_mn Tau(mn,ab)W(mn,ij) +
C    1/2 Sum_ef Tau(ij,ef)W(ab,ef)           (DRLAD)
C
      CALL PCCD_L2INL2(WORK(INEXT),MAXCOR,IUHF,0)

#ifdef _DEBUG_LVL0
       call pccd_check_t2(work(inext),Maxcor,Iuhf)
#endif

C Denominator weigh the T2 increments. 
C
#ifdef _HF_ONLY
      CALL PCCD_SET2NULL(WORK(INEXT),MAXCOR,IUHF,63)
#endif
      IF (PCCD) CALL PCCD_RESET(WORK(INEXT),MAXCOR,IUHF,63)
      CALL PCCD_NEWT2(WORK(INEXT),MAXCOR,IUHF)

#ifdef _DEBUG_LVLM
       call gwinl1(work(Inext),Maxcor,Iuhf)
       call pccd_l2inl1(work(inext),Maxcor,Iuhf)
       call pccd_check_t2(work(inext),Maxcor,Iuhf)
#endif
C
#ifdef _DEBUG_LVL0
      CALL PCCD_DRTSTS(WORK(INEXT),MAXCOR,ICYCLE,IUHF,ICONVG,ICONTL,
     &            .FALSE.,0,'T')
#endif

      IF (IFLAGS(h_IFLAGS_cc_extrapol) .EQ.1)  THEN
         CALL PCCD_DODIIS0(WORK(INEXT),MAXCOR,IUHF,1,ICYCLE,
     &                    ICONVG,ICONTL,.FALSE.,144,61,190,0,90,2,70)
      END IF

      CALL PCCD_CMPENG(WORK(INEXT),MAXCOR,63,ECORR,ENERGY(ICYCLE+1,1))

      IF (ICONVG.NE.0) THEN

         CALL PCCD_DRMOVE(WORK(INEXT),MAXCOR,IUHF,100)
         CALL PCCD_RNABIJ(WORK(INEXT),MAXCOR,IUHF,'L')

      END IF
      IF (ICONVG.EQ.0) THEN

         CALL PCCD_CMPENG(WORK(INEXT),MAXCOR,146,ECORR,ENERGY(1,1))
         CALL PCCD_AMPSUM(WORK(INEXT),MAXCOR,IUHF,0,.FALSE.,'T')

         CALL PCCD_RFINISH(ICYCLE+1,"L")
         CALL PCCD_DDMPTGSS(WORK(INEXT), MAXCOR/IINTFP, IUHF, 0,
     &                      'TGUESS  ')
         WRITE(6,1020)
         RETURN
 1020 FORMAT(/,77('-'),/,32X,'Exiting pccd_ldriver',/,77('-'),/)

      ENDIF 
C
      ICYCLE = ICYCLE +  1
      IF (ICYCLE .LE. NCYCLE) GOTO 100
       
      IF (ICONVG .NE. 0) THEN
         CALL PCCD_AMPSUM(WORK(INEXT),MAXCOR,IUHF,0,.FALSE.,'T')
         Write(6,"(2a)") "The pCCD equations did not converge",
     &                     " in alloted number of cyclces."
      ENDIF 
C
      RETURN
      END
