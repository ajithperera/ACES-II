      SUBROUTINE DRVAHF(SCR, NATOM, NORBS, NUMATOM)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      DOUBLE PRECISION EFG(6)
      CHARACTER*2 LABELS(6)
      DIMENSION SCR(4*NORBS*NORBS), NATOM(NUMATOM)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      DATA LABELS /'XX','YY','ZZ','XY','XZ','YZ'/
      IMXATM=50
      IONE=1
      I000=1
      I010=I000+NORBS*NORBS
      I020=I010+NORBS*NORBS
      I030=I020+NORBS*NORBS
      I040=I030+NORBS*NORBS
      IATOM=0
      IRWND=0
      CALL GETREC(20,'JOBARC','MAP2ZMAT',NUMATOM,NATOM)
30    CALL SEEKLB('   FXX  ',IERR,IRWND)
      IF(IERR.NE.0)RETURN
      IRWND=1
      IF(IATOM.EQ.0)WRITE(LUOUT,1000)
      IATOM=IATOM+1
      CALL COMPPR(EFG(1), SCR(I030), SCR(I010), NORBS, .FALSE.)
      CALL SEEKLB('   FYY  ',IERR,IRWND)
      CALL COMPPR(EFG(2), SCR(I030), SCR(I010), NORBS, .FALSE.)
      CALL SEEKLB('   FZZ  ',IERR,IRWND)
      CALL COMPPR(EFG(3), SCR(I030), SCR(I010), NORBS, .FALSE.)
      CALL SEEKLB('   FXY  ',IERR,IRWND)
      CALL COMPPR(EFG(4), SCR(I030), SCR(I010), NORBS, .FALSE.)
      CALL SEEKLB('   FXZ  ',IERR,IRWND)
      CALL COMPPR(EFG(5), SCR(I030), SCR(I010), NORBS, .FALSE.)
      CALL SEEKLB('   FYZ  ',IERR,IRWND)
      CALL COMPPR(EFG(6), SCR(I030), SCR(I010), NORBS, .FALSE.)

1000  FORMAT(T3,' Dipole-dipole contributions to HFS at',
     &          ' atomic centers ')
      WRITE(LUOUT,100)NATOM(IATOM)
      WRITE(LUOUT,101)(LABELS(I),EFG(I),I=1,6)
101   FORMAT((T6,3(1X,A3,' = ',F15.10)))
100   FORMAT(T30,'Z-matrix center ',I3,':')
      GOTO 30
CJDW  1/8/98. Commented out next line to stop f90 warning.
C     RETURN
      END 