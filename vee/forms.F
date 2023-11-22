
      SUBROUTINE FORMS(LENGTH,NDIM,A,TMP,COEFF,ILIST1,ILIST2,
     &                 IOLDEST,MAXORD)
C
C CALCULATES FULL VECTOR IN TERMS OF A SET OF N BASIS VECTORS
C AND N EXPANSION COEFFICIENTS.  BASIS VECTORS ARE STORED ON
C A LIST.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(LENGTH),TMP(LENGTH),COEFF(NDIM)
C
      IGET(I)=1+MOD(IOLDEST+MAXORD-I,MAXORD+1)
C
      CALL ZERO(A,LENGTH)
      DO 10 I=1,NDIM
       CALL GETLST(TMP,IGET(I),1,1,ILIST1,ILIST2)
#ifdef _DEBUG_LVL0
       Write(6,"(a,i4)") "Checksum of basis vectors in list",ilist2
       call checksum("@-forms          :",tmp,length,s)
       write(6,"(a,3x,F15.7)") "The coefficients            :",coeff(i)
#endif
       CALL SAXPY (LENGTH,COEFF(I),TMP,1,A,1)
10    CONTINUE
C      
      RETURN
      END 
