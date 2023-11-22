      SUBROUTINE SYMMET5(A,B,NUM)
C
C   THIS ROUTINE SYMMETRIZES A GIVEN MATRIX A
C
C     A(PQ) = 1/2 ( A(PQ)+A(QP))
C
C   WHERE A IS A SYMMETRY PACKED MATRIX AND
C   NUM THE CORRESPONDING POPULATION VECTOR
C
C   THE SYMMETRIZED ARRAY IS RETURNED IN A, B IS ONLY USED
C   AS A SCRATCH ARRAY. IT IS POSSIBLE TO SET UP A ROUTINE WHICH
C   DOES NOT USE A SCRATCH ARRAY, BUT VECTORIZATION SHOULD HERE 
C   SLIGHTLY BETTER
C
CEND
C
C CODED JULY/90  JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD
      DIMENSION A(1),B(1),NUM(8)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      DATA HALF /0.5D+0/
C
      IOFF=0
      DO 1000 IRREP=1,NIRREP
C
       NUMI=NUM(IRREP)
C
       DO 100 I=1,NUMI
       DO 100 J=1,NUMI
        B(IOFF+(I-1)*NUMI+J)=HALF*(A(IOFF+(I-1)*NUMI+J)
     &                        +A(IOFF+(J-1)*NUMI+I))
100    CONTINUE
C
       IOFF=IOFF+NUMI*NUMI
1000  CONTINUE
C
C  COPY B BACK TO A
C
      DO 2000 I=1,IOFF 
       A(I)=B(I)
2000  CONTINUE
C
      RETURN
      END