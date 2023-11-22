C
      SUBROUTINE OVLAP(L, M, N, GA, V, NT, D)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C &  Calculate electric field gradient with out Fermi contact   &
C &  contribution. Hacked up version from vprop.f.              &
C &  04/94 Ajith                                                &
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C 
      COMMON /PRTBUG/ IPRTFS  
      COMMON /VARS/ FNU(20),DP(20),FN(15),FD(15),PI,DFTR(11)
      DIMENSION V(10), D(10)
C
      NT = 1
      V(1) = OLAP(L, M, N, GA)
      NT = 3
      RETURN
      END
