C     
      SUBROUTINE GETTYPE (IL, IM, IN, ITYP)
C
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C &  Returns type of the function (eg. S=0, PX=1, PY=2, PZ=3) &
C &  Coded by Ajith 08/93                                     &
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXANG = 7)
      INTEGER INDEX(MAXANG, MAXANG, MAXANG)
C
C
      II = 0
      DO 10 LVAL = 0, (MAXANG -1)
         DO 20 L = LVAL, 0, -1
            LEFT = LVAL - L
            DO 30 M = LEFT, 0, -1
               N = LEFT - M
               II = II + 1
               INDEX(L+1, M+1, N+1) = II
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
C
C
C Add the offset and Handle the case of negative IL, IM, IN.
C see subroutine "OrbPara" for comments about negative
C IL, IM and IN.
C
      IF (IL .LT. 0) THEN 
         IIL = IL + 2
      ELSE
         IIL = IL + 1
      ENDIF
C         
      IF (IM .LT. 0) THEN 
         IIM = IM + 2
      ELSE
         IIM = IM +1
      ENDIF
C      
      IF (IN .LT. 0) THEN
         IIN = IN + 2
      ELSE
         IIN = IN + 1
      ENDIF
C
      ITYP = INDEX(IIL, IIM, IIN)
C     
      RETURN
      END
