    
     
      SUBROUTINE ZERLST(Z,NSIZ,NUMDIS,NCACHE,LISTA,LISTB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Z(NSIZ)
      CALL ZERO(Z,NSIZ)
      CALL PUTLST(Z,1,NUMDIS,NCACHE,LISTA,LISTB)
      RETURN
      END  
