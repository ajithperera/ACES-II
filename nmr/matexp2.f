      SUBROUTINE MATEXP2(NUM,A,B)
C
C     THIS ROUTINE EXPANDS THE A COMPRESSED MATRIX A(P,Q)
C     P >= Q TO AN ARRAY A(PQ) WITH P,Q. NOTE THIS ROUTINE 
C     EXPECTS THAT THE ARRAY A IS SYMMETRY PACKED
C
C     INPUT : NUM ......  POPULATION VECTOR FOR I AND J
C             A     ....  THE MATRIX A
C
C     OUTPUT : B .......  THE EXPANDED MATRIX A
C
CEND
C
C  CODED JG JAN/91
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD
      DIMENSION A(1),B(1)
C
      IND(J,I)=(J*(J-1))/2+I
C 
      DATA ZERO /0.0D0/
C
C     LOOP OVER ORBITALS, BUT ALSO IN BACKWARD ORDER
C
      DO 100 J=1,NUM
       DO 100 I=1,J
        IND1=IND(J,I)
        IND2=(J-1)*NUM+I
        IND3=(I-1)*NUM+J
        B(IND2)=A(IND1)
        B(IND3)=A(IND1)
100   CONTINUE
      RETURN
      END
