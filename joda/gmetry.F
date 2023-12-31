      SUBROUTINE GMETRY(We_havegeom, NOSILENT)
C
C     ROUTINE TO CALCULATE CARTESIAN COORDINATES FROM INTERNAL
C     COORDINATE REPRESENTATION.  SOME OF THIS HAS BEEN LIFTED FROM
C     PRDDO, ALTHOUGH SOME IMPROVEMENTS HAVE BEEN MADE.  UP TO 50
C     ATOMS ALLOWED.  CONNECTIVITY OF FIRST THREE MUST BE 1-2-3 IN
C     INTERNAL COORDINATE REP.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "mxatms.par"
#include "coord.com"
#include "const.h"
      LOGICAL OPTRES,FIRST_FIX
      INTEGER NA(MXATMS),NB(MXATMS),NC(MXATMS)
      DIMENSION A(3,MXATMS),Q0(3*MXATMS)
C
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /FLAGS/ IFLAGS(100),IFLAGS2(500)
C     Main OPTIM control data
C     IPRNT   Print level - not used yet by most routines
C     INR     Step-taking algorithm to use
C     IVEC    Eigenvector to follow (TS search)
C     IDIE    Ignore negative eigenvalues
C     ICURVY  Hessian is in curviliniear coordinates
C     IMXSTP  Maximum step size in millibohr
C     ISTCRT  Controls scaling of step
C     IVIB    Controls vibrational analysis
C     ICONTL  Negative base 10 log of convergence criterion.
C     IRECAL  Tells whether Hessian is recalculated on each cyc
C     INTTYP  Tells which integral program is to be used
C              = 0 Pitzer
C              = 1 VMol
C     XYZTol  Tolerance for comparison of cartesian coordinates
C
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol
C
      COMMON/RESTART2/OPTRES
      COMMON/FIXED/FIRST_FIX 
C
      LOGICAL NOT_IN_BOUND, We_havegeom, NOSILENT
#include "io_units.par"
C
      IONE=1
      CALL GETREC(-1,'JOBARC','PASS1   ',IONE,IX)
CSSS      CALL GETREC(-1,'JOBARC', 'VIB_POPT',IONE, ipost_vib)
      If (IPrnt .ge. 20 .AND. NOSILENT) then
         Write (LuOut, *) 'GMETRY starting with R'
         Write (LuOut,'((I3,1X,f12.6))') (i,R(i),i=1,NX)
      EndIf
C
#ifdef _DEBUG_LVL0
      Print*, "In Gmetry; ncycle, ix, and ipost_vib"
      Print*, ncycle, is, ipost_vib
#endif
      ATOR  = DACOS(-1.D0)/180.D0
      TODEG = 1.0D0/ATOR 
      ATOB  = CONV_a0_a 

      IF (ncycle .NE. 0. or. ix.ne.0) THEN
C
C DECOMPRESS R
C
         If (IPrnt .ge. 20 .AND. NOSILENT) Write (LuOut, '(a)')
     $      '@GMETRY-I, Decompressing R.'
         CALL USQUSH(R,NXM6)
         IF(.NOT.OPTRES .AND. NOSILENT)  WRITE(LuOut,78)
 78      FORMAT(' Updating structure...')
      ELSE
         If (IPrnt .ge. 20 .AND. NOSILENT) Write (LuOut,*)
     $      'GMETRY: Converting to radians & bohr.',ator,atob
         DO 19 IX=8,NX-1,3
            IF(IX.NE.8)R(IX+1)=R(IX+1)*ATOR
 19         R(IX)=R(IX)*ATOR
         IF(IFLAGS(78).EQ.0)THEN
          DO 18 IX=4,NX-2,3
 18          R(IX)=R(IX)/ATOB
         ENDIF
      ENDIF
C
      If (IPrnt .ge. 20 .AND. NOSILENT) then
         Write (LuOut, *) 'GMETRY using R vector'
         Write (LuOut,'((I3,1X,f12.6))') (i,R(i),i=1,NX)
      EndIf

      CALL GEN_CART_COORD(Q0, NOSILENT)
C
C DEAL WITH CASE WHERE ANGLES OR DIHEDRALS NOT WITHIN BOUNDS HERE.
C
      CALL ZERO(Q0,3*MXATMS)
      CALL DCOPY(NX,R,1,Q0,1)
      CALL XTOR(R,0)
      CALL USQUSH(R,NXM6)
      NOT_IN_BOUND = .FALSE.
      DO 107 I=1,NX
       IF(DABS(DABS(R(I))-DABS(Q0(I))).GT.1.D-4)THEN
        IF(NOSILENT) WRITE(6,8331)I,Q0(I)*TODEG,R(I)*TODEG 
 8331   FORMAT('@GMETRY-W, Internal coordinate #',i3,
     &' not within bounds:',
     &/,t3,' Value was: ',f10.5,' and has been changed to ',f10.5,'.')
        NOT_IN_BOUND = .TRUE.
        CALL GETREC(0,"JOBARC","DIH_FIX ",LENGTH,JUNK)
        IF (LENGTH .LE. 0) CALL FIXDIH(I,NX) 
       ENDIF
 107  CONTINUE
C
C Since their inception, the previous 12 lines, which were meant to
C correct poorly constructed internal coordinates, were doing more harm
C than good for geometry optimizations. I assume that this was left
C unnoticed for so many years simply because every one wrote decent
C Z-matrices. Then come ZMAT files generated automatically by MOPAC (or
C HyperChem?). Those matrices go through the 12 lines and trigger the
C "not within bound" flag, which causes joda to use its own internally-
C chosen values. The problem is that the Cartesian coordinates still
C correspond to the "bad" internal coordinates. For single point
C calculations, this is irrelevant since the internal coordinates have
C no role beyond generating the Cartesian coordinates (so why do we even
C bother for SP calcs?). However, during an optimization, we go back
C and forth between Cartesians and internals. As a result, when the
C condition that ((DABS(DABS(R(I))-DABS(Q0(I))).GT.1.D-4) is satisfied,
C the gradients are calculated for Cartesians that do NOT correspond to
C the internals that are supposed to be updated. This is a serious error
C and might partially explain why some TS searches were "meandering". To
C fix this, I created a new subroutine called GEN_CART_COORD which takes
C internal coordinates and connectivities (from ACES II common blocks)
C and generates ACES II Cartesians (nothing else). Ajith Perera, 04/2005
C
C Only one out-of-bound is all it takes to regenerate Cartesians.
      
      IF (NOT_IN_BOUND)  THEN
         CALL GEN_CART_COORD(Q0, NOSILENT)
      ENDIF 
C
#ifdef _DEBUG_LVLM
      Call geomout
#endif
      RETURN
      END
