      SUBROUTINE DRVEFG(SCR, NATOM, NORBS, NUMATOM)
c
c Slightly modified by M. Nooijen to print out maximum principal value
c and asymmetry parameter
c
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION EFG(6), efg_tensor(3,3), evec(3,3),
     $     eval(3), iloc(3), efg_max, efg_asym
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
      CALL COMPPR(EFG(1), SCR, SCR(I030), NORBS, .TRUE.)
      efg_tensor(1,1) = EFG(1)
#ifdef _DEBUG_LVLM
      Write(6,"(a)") " The EFG-XX " 
      call output(SCR(I030), 1, norbs, 1, norbs, norbs, norbs, 1)
#endif
      CALL SEEKLB('   FYY  ',IERR,IRWND)
      CALL COMPPR(EFG(2), SCR, SCR(I030), NORBS, .TRUE.)
      efg_tensor(2,2) = EFG(2)
#ifdef _DEBUG_LVLM
      Write(6,"(a)") " The EFG-YY " 
      call output(SCR(I030), 1, norbs, 1, norbs, norbs, norbs, 1)
#endif
      CALL SEEKLB('   FZZ  ',IERR,IRWND)
      CALL COMPPR(EFG(3), SCR, SCR(I030), NORBS, .TRUE.)
      efg_tensor(3,3) = EFG(3)
#ifdef _DEBUG_LVLM
      Write(6,"(a)") " The EFG-ZZ " 
      call output(SCR(I030), 1, norbs, 1, norbs, norbs, norbs, 1)
#endif
      CALL SEEKLB('   FXY  ',IERR,IRWND)
      CALL COMPPR(EFG(4), SCR, SCR(I030), NORBS, .TRUE.)
      efg_tensor(1,2) = EFG(4)
      efg_tensor(2,1) = EFG(4)
#ifdef _DEBUG_LVLM
      Write(6,"(a)") " The EFG-XY " 
      call output(SCR(I030), 1, norbs, 1, norbs, norbs, norbs, 1)
#endif
      CALL SEEKLB('   FXZ  ',IERR,IRWND)
      CALL COMPPR(EFG(5), SCR, SCR(I030), NORBS, .TRUE.)
      efg_tensor(1,3) = EFG(5)
      efg_tensor(3,1) = EFG(5)
#ifdef _DEBUG_LVLM
      Write(6,"(a)") " The EFG-XZ " 
      call output(SCR(I030), 1, norbs, 1, norbs, norbs, norbs, 1)
#endif
      CALL SEEKLB('   FYZ  ',IERR,IRWND)
      CALL COMPPR(EFG(6), SCR, SCR(I030), NORBS,. TRUE.)
      efg_tensor(3,2) = EFG(6)
      efg_tensor(2,3) = EFG(6)
#ifdef _DEBUG_LVLM
      Write(6,"(a)") " The EFG-YZ " 
      call output(SCR(I030), 1, norbs, 1, norbs, norbs, norbs, 1)
#endif
      call process_efg_tensor(efg_tensor, evec, eval, efg_max,
     $     efg_asym)

1000  FORMAT(T3,' Electric field gradient at atomic centers ')
      WRITE(LUOUT,100)NATOM(IATOM)
      WRITE(LUOUT,101)(LABELS(I),EFG(I),I=1,6)
      WRITE(LUOUT,102), Natom(iatom), efg_max, efg_asym
101   FORMAT((T6,3(1X,A3,' = ',F15.10)))
100   FORMAT(T30,'Z-matrix center ',I3,':')
 102  FORMAT(T4, 'Principal axis EFG: ', ' atom = ', I3,
     $     ' Max. = ', F15.10, ' Asym. = ', F15.10)
      GOTO 30
C     RETURN
      END 
