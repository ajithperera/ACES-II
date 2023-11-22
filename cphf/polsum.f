
      SUBROUTINE POLSUM(SCF)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SCF
C
      COMMON /HYPER/BETA(3,3,3)
      COMMON /POLAR/  POLARS(3,3), POLFLT(3,3)
C
      IF(SCF) THEN
       CALL HEADER('Static polarizability',-1)
      ELSE
       CALL HEADER('SCF static polarizability',-1)
      ENDIF 
      CALL POLPRI(POLARS,'au')
C
      IF(SCF) THEN
       CALL HEADER('Static hyperpolarizability',-1)
      ELSE
       CALL HEADER('SCF static hyperpolarizability',-1)
      ENDIF
      CALL HYPPRI(BETA,'au')
C
      RETURN
      END
