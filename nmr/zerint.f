      SUBROUTINE ZERINT(ICORE,MAXCOR)  
      IMPLICIT INTEGER(A-Z)
      DIMENSION ICORE(MAXCOR)
C
      CALL ZEROLIST(ICORE,MAXCOR,316)
      CALL ZEROLIST(ICORE,MAXCOR,330)
      CALL ZEROLIST(ICORE,MAXCOR,325)
      CALL ZEROLIST(ICORE,MAXCOR,321)
      CALL ZEROLIST(ICORE,MAXCOR,310)
C
      RETURN
      END
