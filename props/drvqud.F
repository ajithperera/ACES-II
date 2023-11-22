C
      SUBROUTINE DRVQUD(SCR,MAXCOR,NORBS)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "flags.h"
C
      DIMENSION SCR(MAXCOR)
      DOUBLE PRECISION QMOM(6), XPROP(19)
      CHARACTER*2 LABELS(6)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS2/ IFLAGS2(500)
      DATA LABELS /'XX','YY','ZZ','XY','XZ','YZ'/
      I010=1+NORBS*NORBS*IINTFP
      CALL SEEKLB('   QXX  ',IERR,0)
      IF(IERR.NE.0)RETURN
      CALL COMPPR(QMOM(1),SCR,SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('   QYY  ',IERR,1)
      CALL COMPPR(QMOM(2),SCR,SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('   QZZ  ',IERR,1)
      CALL COMPPR(QMOM(3),SCR,SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('   QXY  ',IERR,1)
      CALL COMPPR(QMOM(4),SCR,SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('   QXZ  ',IERR,1)
      CALL COMPPR(QMOM(5),SCR,SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('   QYZ  ',IERR,1)
      CALL COMPPR(QMOM(6),SCR,SCR(I010),NORBS,.TRUE.)
      WRITE(LUOUT,100)
100   FORMAT(T3,' Components of electric quadrupole moment ')
      WRITE(LUOUT,101)(LABELS(I),QMOM(I),I=1,6)
101   FORMAT((T6,3(1X,A3,' = ',F15.10)))
c
c write properties to jobarc
c
      if (.true.) then
      if (iflags2(h_IFLAGS2_ffield) .ne. 0) then
      call getrec(20,'JOBARC','MULTMOM ',
     $              19*iintfp, xprop)
      do i = 1, 6
         xprop(i+3) = qmom(i)
      enddo
      call putrec(20,'JOBARC','MULTMOM ',
     $              19*iintfp, xprop)
      endif
      endif
c
      RETURN
      END 
