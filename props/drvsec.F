      SUBROUTINE DRVSEC(SCR,MAXCOR,NORBS)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "flags.h"
C
      DIMENSION SCR(MAXCOR)
      DOUBLE PRECISION SECMOM(6), XPROP(19)
      CHARACTER*2 LABELS(6)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS2/ IFLAGS2(500)
      DATA LABELS /'XX','YY','ZZ','XY','XZ','YZ'/
      I010=1+NORBS*NORBS
      CALL SEEKLB('   2XX  ',IERR,0)
      IF(IERR.NE.0)RETURN
      CALL COMPPR(SECMOM(1), SCR, SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('   2YY  ',IERR,1)
      CALL COMPPR(SECMOM(2), SCR, SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('   2ZZ  ',IERR,1)
      CALL COMPPR(SECMOM(3), SCR, SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('   2XY  ',IERR,1)
      CALL COMPPR(SECMOM(4), SCR, SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('   2XZ  ',IERR,1)
      CALL COMPPR(SECMOM(5), SCR, SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('   2YZ  ',IERR,1)
      CALL COMPPR(SECMOM(6), SCR, SCR(I010),NORBS,.TRUE.)
      WRITE(LUOUT,100)
100   FORMAT(T3,' Components of second moment ')
      WRITE(LUOUT,101)(LABELS(I),SECMOM(I),I=1,6)
101   FORMAT((T6,3(1X,A3,' = ',F15.10)))
c
c write properties to jobarc
c
      if (.false.) then
c These are not the correct second order multipoles. See drvqud.
      if (iflags2(h_IFLAGS2_ffield) .ne. 0) then
      call getrec(20,'JOBARC','MULTMOM ',
     $              19*iintfp, xprop)
      do i = 1, 6
         xprop(i+3) = secmom(i)
      enddo
      call putrec(20,'JOBARC','MULTMOM ',
     $              19*iintfp, xprop)
      endif
      endif
c
      RETURN
      END 
