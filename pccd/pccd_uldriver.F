#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */
      SUBROUTINE PCCD_ULDRIVER(ICORE,MAXCOR,IUHF)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION ICORE(MAXCOR)
      DIMENSION ECORR(3)
      INTEGER RELCYC

      LOGICAL NONHF
      LOGICAL CIS,EOM
      LOGICAL DO_TAU
      LOGICAL PCCD,CCD,LCCD
C
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

      I0      = 1 
      INEXT   = 1
      RLECYC  = 0
      IREF    = 1
      IBWPASS = 0
      IRREPX  = 1
      ICONTL  = IFLAGS(h_IFLAGS_cc_conv)
      NCYCLE  = IFLAGS(h_IFLAGS_cc_maxcyc)
      NKEEP   = IFLAGS(h_IFLAGS_cc_exporder)
      IEVERY  = IFLAGS(h_IFLAGS_tamp_sum)
      IPRTENG = 1

C-----------------------------------------------------------------------

      CALL PCCD_INILAM(IUHF)
      CALL ZERO(ENERGY,1000)

      IF (LCCD) THEN
         Write(6,"(a)") "      Lambda is constructed for LCCD."
         CALL MVSYM(ICORE(I0),MAXCOR,IUHF,IRREPX,46,146)
         CALL GETLST(ICORE(I0),1,1,1,1,90)
         CALL PUTLST(ICORE(I0),1,1,1,1,190)
         CALL LGENINT(ICORE(I0),MAXCOR,IUHF)
         WRITE(6,1020)
         RETURN
      ENDIF

      WRITE(6,801)
      CALL AMPSUM(ICORE(i0),MAXCOR,IUHF,100,SING,'L')
801   FORMAT('  Initial lambda amplitudes: ')
      WRITE(6,800)
  800 FORMAT(T3,' Correlation energies computed from initial ',
     &          'L amplitudes: ')
      CALL CMPENG2(ICORE(i0),MAXCOR,143,0,ECORR,ENERGY(1,1),
     &              ENERGY(1,2),IUHF,IPRTENG)
C
      IF (IFLAGS(h_IFLAGS_cc_extrapol) .EQ.1) THEN
         CALL PCCD_DIISLST(IRREPX,IUHF,.FALSE.,"L")
      ENDIF 
C
      CALL RNABIJ(ICORE(i0),MAXCOR,IUHF,'L')
      CALL STLLST(ICORE(i0),MAXCOR,IUHF)
C
C CALCULATE FOR CCSD THE W(AB,EF) INTERMEDIATES
C
C
      CALL FORMWL(ICORE(i0),MAXCOR,IUHF) 
      CALL WTWTW(ICORE(i0),MAXCOR,IUHF)

      TERM1 = .TRUE.
      TERM2 = .TRUE.
      TERM3 = .TRUE.
      TERM4 = .TRUE.
      TERM5 = .TRUE.
      TERM6 = .TRUE.

      DO_TAU = .TRUE.
      CALL FORMW4(ICORE(i0),MAXCOR,IUHF,TERM1,TERM2,TERM3,
     &            TERM4,TERM5,TERM6,DO_TAU,.TRUE.)

      TERM1=.TRUE.
      TERM2=.TRUE.
      TERM3=.TRUE.
      TERM4=.TRUE.
      TERM5=.TRUE.
      TERM6=.TRUE.
      CALL FORMW5(ICORE(i0),MAXCOR,IUHF,TERM1,TERM2,TERM3,
     &            TERM4,TERM5,TERM6,DO_TAU,.TRUE.)

c************************************************************
c from now on lists 33-39 (newly created taus) are on 233-239
c*************************************************************

      ICYCLE = 1
 100  CONTINUE

      RLECYC=RLECYC+1
C
C FILL LAMBDA(2) INCREMENTS WITH <IJ||AB> INTEGRALS
C
C AND ZERO THE LAMBDA(1) INCREMENTS (ONLY CCSD AND QCISD)
C
      CALL INITIN(ICORE(i0),MAXCOR,IUHF)
CSSS      CALL INITSN(ICORE(i0),MAXCOR,IUHF)
C
C GENERATE V AND G INTERMEDIATES.
C
#ifdef _DEBUG_LVL0
      call check_leom(icore(i0),Maxcor,Iuhf)
#endif 

      CALL LGENINT(ICORE(i0),MAXCOR,IUHF)

#ifdef _DEBUG_LVL0
      call check_leom(icore(i0),Maxcor,Iuhf)
#endif 
C
C COMPUTE F AND G INTERMEDIATE CONTRIBUTION TO LAMBDA(2) INCREMENT.
C
      CALL F1INL2(ICORE(i0),MAXCOR,IUHF)
      CALL F2INL2(ICORE(i0),MAXCOR,IUHF)

#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals before GinL2"
      call check_leom(icore(i0),Maxcor,Iuhf)
#endif 

      CALL PCCD_G1INL2U(ICORE(i0),MAXCOR,IUHF)
      CALL PCCD_G2INL2U(ICORE(i0),MAXCOR,IUHF)

#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals after GinL2"
      call check_leom(icore(i0),Maxcor,Iuhf)
#endif
C
C DO W INTERMEDIATE CONTRIBUTION TO T2 EQUATION (SAME AS THIRD ORDER CODE).
C
      CALL PCCD_L2INL2U(ICORE(i0),MAXCOR,IUHF,0)
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals after l2lad"
      call check_leom(icore(i0),Maxcor,Iuhf)
#endif
      IF (pCCD)  THEN
C Ajith Perera, 10/2021.
         Write(6,"(20x,a)") " ----------Warning------------"
         Write(6,"(3a)")  "  A pCCD or pCCD like calculation is being",
     +                   " perforemed and the off-diagonal"
         write(6,"(a)")   " blocks of T2 is set to zero"
         Write(6,"(20x,a)") " -----------------------------"
         call pccd_reset_vcc(Icore(i0),Maxcor,Iuhf,63)
      ENDIF 
#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals after l2inl2"
      call check_leom(icore(i0),Maxcor,Iuhf)
#endif
C
C DENOMINATOR WEIGHT T2 INCREMENTS TO FORM NEW T2.
C
      CALL NEWT2(ICORE(i0),MAXCOR,IUHF)

#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals after newT2"
      call check_leom(icore(i0),Maxcor,Iuhf)
#endif
#ifdef _DEBUG_LVL0
      CALL DRTSTS(ICORE(i0),MAXCOR,ICYCLE,IUHF,ICONVG,ICONTL,
     &              SING,100,'L')
#endif 
      CALL CMPENG2(ICORE(i0),MAXCOR,60,2,ECORR,ENERGY(ICYCLE+1,1),
     &               ENERGY(ICYCLE+1,2),IUHF,IPRTENG)

      IF (IFLAGS(h_IFLAGS_cc_extrapol) .EQ.1) THEN
        CALL PCCD_DODIIS0(ICORE(i0),MAXCOR/IINTFP,IUHF,1,ICYCLE,
     &                    ICONVG,ICONTL,.FALSE.,144,61,190,0,90,
     &                    2,70)
      ENDIF

      IF (ICONVG .NE. 0)  THEN
         CALL DRMOVE(ICORE(i0),MAXCOR,IUHF,100,SING)
         CALL RNABIJ(ICORE(i0),MAXCOR,IUHF,'L')
      ENDIF 

      IF(ICONVG.EQ.0) THEN

         CALL CMPENG2(ICORE(i0),MAXCOR,60,2,ECORR,ENERGY(ICYCLE+1,1),
     &                ENERGY(ICYCLE+1,2),IUHF,IPRTENG)
         CALL AMPSUM(ICORE(i0),MAXCOR,IUHF,0,.FALSE.,'L')

CSSS         CALL FINISH(ICYCLE+1,PCCD,.FALSE.,.FALSE.)
         CALL PCCD_UFINISH(ICYCLE+1,"L")
         CALL PCCD_DDMPTGSS(ICORE(i0), MAXCOR/IINTFP, IUHF, 0,
     &                      'TGUESS  ')
         WRITE(6,1020)
 1020    FORMAT(/,77('-'),/,32X,'Exiting pccd_ldriver',/,77('-'),/)
         RETURN
      ENDIF

      ICYCLE = ICYCLE +  1
      IF (ICYCLE .LE. NCYCLE) GOTO 100

      IF (ICONVG .NE. 0) THEN
         CALL AMPSUM(ICORE(i0),MAXCOR,IUHF,0,.FALSE.,'T')
         Write(6,"(2a)") "The pCCD equations did not converge",
     &                     " in alloted number of cyclces."
      ENDIF

      RETURN
      END
