C
      SUBROUTINE DUMPJ(VLIST, NOC)
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C & Write atom ordering in vprop.f for future use in NMR     &
C & coupling constant calculations.                          &
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
C Here there is a hard limit on number of centers (100).
C Did not try to fix these hard limits since this program
C going to be obsolete.
C
      DIMENSION VLIST(NOC, 10), ITMP(100), CENTER(100, 3)
      DO 10 I= 1, NOC
         ITMP(I) = VLIST(I, 4) 
 10   CONTINUE
C
      DO 20 II = 1, NOC
         DO 30 JJ = 1, 3
            CENTER(II, JJ) = VLIST(II, JJ)
 30      CONTINUE
 20   CONTINUE
C         
      CALL PUTREC (20, 'JOBARC', 'VMOLORDR', NOC, ITMP)
      CALL PUTREC (20, 'JOBARC', 'VMORDCOR', 100*3*IINTFP, CENTER)
C
      RETURN
      END      
