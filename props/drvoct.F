
      SUBROUTINE DRVOCT(SCR,MAXCOR,NORBS)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
#include "flags.h"
C
      DIMENSION SCR(MAXCOR)
      DOUBLE PRECISION OCTMOM(10), XPROP(19)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS2/ IFLAGS2(500)
      CHARACTER*3 LABELS(10)
      DATA LABELS /'XXX','YYY','ZZZ','XXY','XXZ','XYY','YYZ',
     &             'XZZ','YZZ','XYZ'/
      I010=1+NORBS*NORBS
      CALL SEEKLB('  OXXX  ',IERR,0)
      IF(IERR.NE.0)RETURN
      CALL COMPPR(OCTMOM(1), SCR, SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('  OYYY  ',IERR,1)
      CALL COMPPR(OCTMOM(2), SCR, SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('  OZZZ  ',IERR,1)
      CALL COMPPR(OCTMOM(3), SCR, SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('  OXXY  ',IERR,1)
      CALL COMPPR(OCTMOM(4), SCR, SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('  OXXZ  ',IERR,1)
      CALL COMPPR(OCTMOM(5), SCR, SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('  OXYY  ',IERR,1)
      CALL COMPPR(OCTMOM(6), SCR, SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('  OYYZ  ',IERR,1)
      CALL COMPPR(OCTMOM(7), SCR, SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('  OXZZ  ',IERR,1)
      CALL COMPPR(OCTMOM(8), SCR, SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('  OYZZ  ',IERR,1)
      CALL COMPPR(OCTMOM(9), SCR, SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('  OXYZ  ',IERR,1)
      CALL COMPPR(OCTMOM(10), SCR, SCR(I010),NORBS,.TRUE.)
      WRITE(LUOUT,100)
100   FORMAT(T3,' Components of electric octopole moment ')
      WRITE(LUOUT,101)(LABELS(I),OCTMOM(I),I=1,10)
101   FORMAT((T6,3(1X,A3,' = ',F15.10)))
c
c
c write properties to jobarc
c
      if (iflags2(h_IFLAGS2_ffield) .ne. 0) then
      call getrec(20,'JOBARC','MULTMOM ',
     $              19*iintfp, xprop)
      do i = 1, 10
         xprop(9+i) = octmom(i)
      enddo
      call putrec(20,'JOBARC','MULTMOM ',
     $              19*iintfp, xprop)
      endif
c
      RETURN
      END 
