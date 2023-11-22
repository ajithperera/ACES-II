      SUBROUTINE EXPND2A(WPACK,WFULL,NDIM)
C
C  This routine expands a triangularly packed vector of numbers
C  into an antisymmetric square matrix.
C
C   WPACK((NDIM*(NDIM+1))/2) ==> WFULL(NDIM,NDIM)
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WPACK((NDIM*(NDIM+1))/2),WFULL(NDIM,NDIM)
      ITHRU=0
      DO 10 I=1,NDIM
CDIR$ IVDEP
       DO 20 J=1,I
        ITHRU=ITHRU+1 
        WFULL(I,J)=WPACK(ITHRU)
        WFULL(J,I)=-WPACK(ITHRU)
20     CONTINUE
10    CONTINUE
      RETURN
      END
