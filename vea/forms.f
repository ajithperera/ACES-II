C
C ***************************************************************
C  SUBROUTINES THAT CALCULATE THINGS FROM ITERATION VECTORS
C  AND DETERMINE NEW ITERATION VECTORS
C ***************************************************************
C
      SUBROUTINE FORMS(LENGTH,NDIM,A,TMP,COEFF,ISIDE,ILIST,
     &                 IOLDEST,MAXORD)
C
C CALCULATES FULL VECTOR IN TERMS OF A SET OF N BASIS VECTORS
C AND N EXPANSION COEFFICIENTS.  BASIS VECTORS ARE STORED ON
C LIST (ISIDE, ILIST)
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(LENGTH),TMP(LENGTH),COEFF(NDIM)
C
      IGET(I)=1+MOD(IOLDEST+MAXORD-I,MAXORD+1)
C
      CALL ZERO(A,LENGTH)
      DO 10 I=1,NDIM
       CALL GETLST(TMP,IGET(I),1,1,ISIDE,ILIST)
       CALL SAXPY (LENGTH,COEFF(I),TMP,1,A,1)
10    CONTINUE
C      
      RETURN
      END 
