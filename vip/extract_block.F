      SUBROUTINE EXTRACT_BLOCK(W2T2FULL,W2T2SUB,NROW,NCOL,NBLOCK,IBEGN)
     +                          

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION W2T2FULL(NROW,NCOL), W2T2SUB(NROW,NBLOCK)

      DATA IONE /1/

      IOFF = IONE 
      DO I = IBEGN, NCOL
         CALL DCOPY(NROW,W2T2FULL(1,I),1,W2T2SUB(1,IOFF),1)
         IOFF = IOFF + 1
      ENDDO 

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "W2/t2 full" 
      Call output(W2t2full,1,nrow,1,ncol,nrow,ncol,1)
      Write(6,"(a)") "W2/t2 sub" 
      Call output(W2t2sub,1,nrow,1,nblock,nrow,nblock,1)
#endif 

      RETURN 
      END
