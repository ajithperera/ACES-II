      SUBROUTINE QU0SUM(QUADMN,QUADME,QUAD0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QUADMN(6),QUADME(6),QUAD0(6)
C
      COMMON/PERT/NTPERT,NPERT(8),IPERT(8),IXYZ(3),
     &            ISCHROTT(7)
      DATA AZERO /0.D0/
C
      IND=0
      DO 100 I=1,3
      DO 100 J=I,3
      IND=IND+1
      IF (IXYZ(I).EQ.IXYZ(J)) THEN
       QUAD0(IND)=QUADMN(IND)-QUADME(IND)
      ELSE
       QUAD0(IND)=AZERO
      END IF
  100 CONTINUE  
C
cmn      CALL HEADER('Total quadrupole moment',-1) 
cmn      CALL QU0PRI(QUAD0) 
C
      RETURN
      END
