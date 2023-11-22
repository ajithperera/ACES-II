C
 
 
      SUBROUTINE UPDVC(CORE,ENERGY,ICYCLE,MAXCYC,IUHF)
C
C UPDATES THE "OLD ENERGY" VECTOR
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CORE(MAXCYC),ENERGY(3)
      CALL GETLST(CORE,1,1,1,9,50)
      CORE(ICYCLE)=ENERGY(3)
      CALL PUTLST(CORE,1,1,1,9,50)
      IF(IUHF.NE.0)THEN
       DO 10 ISPIN=1,2
        CALL GETLST(CORE,1,1,1,ISPIN+6,50)
        CORE(ICYCLE)=ENERGY(ISPIN)
        CALL PUTLST(CORE,1,1,1,ISPIN+6,50)
10     CONTINUE
      ENDIF
      RETURN
      END