C
      SUBROUTINE COMPUTE_T1_4T2(ICORE,MAXCOR,IUHF,DIFF)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION ICORE(MAXCOR)
      LOGICAL DO_HBAR_4LCCSD
C
C This routine computes the full T1 contribution regardless of 
C The nature of T2 (can be CCSD, CCD, MBPT(2),....etc.)
C
#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */
#include "machsp.com"
#include "syminf.com" /* for nirrep */
#include "sym.com"
C
C Obtained the current T1 vector and save it in memory. Note
C at the first-iteration, T1 vector is zero. 
C
#ifdef _DEBUG_LVL0
      CALL AMPSUM(ICORE(1),MAXCOR,IUHF,0,.TRUE.,'T')
#endif 
C
      DO_HBAR_4LCCSD = .FALSE. 
      IT1OLD_AA = 1
      IT1NEW_AA = IT1OLD_AA + (NT(1)+IUHF*NT(2))*IINTFP
      INEXT     = IT1NEW_AA + (NT(1)+IUHF*NT(2))*IINTFP
      MAXCOR    = MAXCOR - INEXT 
C
      CALL GETLST(ICORE(IT1OLD_AA),1,1,1,1,90)
      IT1OLD_BB  = IT1OLD_AA + NT(1)*IINTFP
      IF (IUHF .NE. 0) CALL GETLST(ICORE(IT1OLD_BB),1,1,1,2,90)
C
C Generate W(mbej), W(mnij) W(mbej, F(ea), F(ij) and F(me)
C intermediates.
C
      CALL GENINT(ICORE(INEXT),MAXCOR,IUHF,2,GEN_HBAR_4LCCSD)
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
      CALL FMECONT(ICORE(INEXT),MAXCOR,IUHF,1)
      IF(IUHF.NE.0) CALL FMECONT(ICORE(INEXT),MAXCOR,IUHF,2)
C
C Evaluate -Sum_nf T(n,f)<na||if> contribution to T1.
C
      CALL T1INT1(ICORE(INEXT),MAXCOR,IUHF,1)
      IF(IUHF.NE.0) CALL T1INT1(ICORE(INEXT),MAXCOR,IUHF,2)
C
C Evaluate the following three contributions to T1:
C
C   -1/2 Sum_mef T(im,ef)<ma||ef>  (T2T1AA1, T2T1AB1)
C   -1/2 Sum_men T(nm,ei)<nm||ei>  (T2T1AA2, T2T1AB2)
C
      CALL E3S(ICORE(INEXT),MAXCOR,IUHF,1)
C 
C Denominator weigh the T1 increments. Need to be carefull 
C about settings of the last argument. If it is postive
C the contribution in list 93 3 and 4 (f(a,i)/Dai terms) are
C added. 
C
      CALL NEWT1(ICORE(INEXT),MAXCOR,IUHF,1)
C
C Copy the T1 from 90(3,4) to 90(1,2).
C
      DO ISPIN=1,IUHF+1
          CALL GETLST(ICORE(INEXT),1,1,1,ISPIN+2,90)
          CALL PUTLST(ICORE(INEXT),1,1,1,ISPIN,90)
      ENDDO
C
C At this point, we should have all the t1 contributions 
C computed. So, we can look at the difference between
C the current and previouse T1 vector. 
C 
#ifdef _DEBUG_LVL0
      CALL AMPSUM(ICORE(INEXT),MAXCOR,IUHF,0,.TRUE.,'T')
#endif 
      CALL GETLST(ICORE(IT1NEW_AA),1,1,1,1,90)
      IT1NEW_BB  = IT1NEW_AA + NT(1)*IINTFP
      IF (IUHF .NE. 0) CALL GETLST(ICORE(IT1NEW_BB),1,1,1,2,90)

      LENGTH = (NT(1)+IUHF*NT(2))*IINTFP

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The previous T1 amplitudes"
      call output(ICORE(IT1OLD_AA),1,nt(1),1,1,nt(1),1,1)
      Write(6,"(a)") "The current T1 amplitudes"
      call output(ICORE(IT1NEW_AA),1,nt(1),1,1,nt(1),1,1)
#endif

      CALL SAXPY(LENGTH,-1.0D0,ICORE(IT1OLD_AA),1,ICORE(IT1NEW_AA),1)
      DIFF = FNDLRGAB(ICORE(IT1NEW_AA),LENGTH)

#ifdef _DEBUG_LVL0
      Write(6,"(a,F15.10)") "The largest T1 diff", diff
#endif
     
      RETURN
      END
