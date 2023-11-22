C
      SUBROUTINE COMPUTE_T1_4T2(ICORE,MAXCOR,IUHF,DIFF)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION ICORE(MAXCOR)
      LOGICAL DO_HBAR_4LCCSD,ORDER2,NONHF
C
C This routine computes the full T1 contribution regardless of 
C The nature of T2 (can be CCSD, CCD, MBPT(2),....etc.)
C
      COMMON /NHFREF/  NONHF
C
#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */
#include "machsp.com"
#include "syminf.com" /* for nirrep */
#include "sym.com"
#include "flags.com"
#include "flags2.com"
C
C Obtained the current T1 vector and save it in memory. Note
C at the first-iteration, T1 vector is zero. 
C
#ifdef _DEBUG_LVLM
      CALL AMPSUM(ICORE(1),MAXCOR,IUHF,0,.TRUE.,'T')
#endif 
C
      DO_HBAR_4LCCSD = .FALSE. 
      ORDER2    =  (Iflags(h_IFLAGS_calclevel) .EQ. 1)
      IT1OLD_AA = 1
      IT1NEW_AA = IT1OLD_AA + (NT(1)+IUHF*NT(2))*IINTFP
      INEXT     = IT1NEW_AA + (NT(1)+IUHF*NT(2))*IINTFP
      MAXCOR    = MAXCOR - INEXT 
      LENGTH    = (NT(1)+IUHF*NT(2))*IINTFP
C
      CALL GETLST(ICORE(IT1OLD_AA),1,1,1,1,90)
      IT1OLD_BB  = IT1OLD_AA + NT(2)*IINTFP
      IF (IUHF .NE. 0) CALL GETLST(ICORE(IT1OLD_BB),1,1,1,2,90)

      If (ORDER2) THEN

C Form f(a,e) contribution to F(a,e)
C
         IF (NONHF) THEN
C
            CALL GETLST(ICORE(INEXT),1,1,1,3,92)
            CALL PUTLST(ICORE(INEXT),1,1,1,1,92)
            IF (IUHF .NE. 0) THEN
               CALL GETLST(ICORE(INEXT),1,1,1,4,92)
               CALL PUTLST(ICORE(INEXT),1,1,1,2,92)
            ENDIF 
C
C Form F(a,e) contribution to T1
C
            CALL FEACONT1(ICORE(INEXT),MAXCOR,IUHF)

C Form f(m,i) contribution to F(M,I)

            CALL GETLST(ICORE(INEXT),1,1,1,3,91)
            CALL PUTLST(ICORE(INEXT),1,1,1,1,91)

            IF (IUHF .NE. 0) THEN
               CALL GETLST(ICORE(INEXT),1,1,1,4,91)
               CALL PUTLST(ICORE(INEXT),1,1,1,2,91)
           ENDIF 
C
C Form F(m,i) contribution to T1
C
            CALL FMICONT1(ICORE(INEXT),MAXCOR,IUHF)
C
C form f(m,e) contribution to F(m,e)
C 
            CALL GETLST(ICORE(INEXT),1,1,1,3,93)
            FAI=DDOT(LENGTH,ICORE(INEXT),1,ICORE(INEXT),1)
            IF (FAI .LT. 1.0D-6) THEN
               Write(6, "(a,F10.7,a)") "The norm of f(a,i)", FAI,
     +         " is too small for have 2nd-order T1 contributions."
               Call ERREX
            ENDIF 
            CALL PUTLST(ICORE(INEXT),1,1,1,1,93)
            IF (IUHF .NE. 0) THEN
               CALL GETLST(ICORE(INEXT),1,1,1,4,93)
               CALL PUTLST(ICORE(INEXT),1,1,1,2,93)
           ENDIF 
C
C Evalaute Sum_e T(im,ae)F(a,i) contribution to T1.
C  
           CALL FMECONT(ICORE(INEXT),MAXCOR,IUHF,1)
           IF(IUHF.NE.0) CALL FMECONT(ICORE(INEXT),MAXCOR,IUHF,2)

          ELSE
            Write(6,"(a)") " COPT-MBPT(2) require NON-HF orbitals."
            Call Errex 
          ENDIF 
      Else
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

      Endif 
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
      IT1NEW_BB  = IT1NEW_AA + NT(2)*IINTFP
      IF (IUHF .NE. 0) CALL GETLST(ICORE(IT1NEW_BB),1,1,1,2,90)

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The previous T1 amplitudes"
      call output(ICORE(IT1OLD_AA),1,nt(1),1,1,nt(1),1,1)
      Write(6,"(a)") "The current T1 amplitudes"
      call output(ICORE(IT1NEW_AA),1,nt(1),1,1,nt(1),1,1)
#endif

      CALL SAXPY(LENGTH,-1.0D0,ICORE(IT1OLD_AA),1,ICORE(IT1NEW_AA),1)
      DIFF = FNDLRGAB(ICORE(IT1NEW_AA),LENGTH)

#ifdef _DEBUG_LVL0
      Write(6,"(2x,a,F15.10)") "The largest T1 diff", diff
#endif
     
      RETURN
      END
