      SUBROUTINE WHENFGT(N,ARRAY,INC,TARGET,INDEX,NVAL)
C
C     Returns the locations of all elements in a real ARRAY
C          that are greater than a real TARGET.
C     Martin J. McBride.  7/17/85.
C     General Electric CRD, Information System Operation.
C
      INTEGER N,INC,INDEX(1),NVAL,IX,I
      DOUBLE PRECISION ARRAY(1),TARGET

      IF (N .LT. 1) RETURN
      IX = 1
      NVAL = 0
      IF (INC .LT. 0) IX = (-INC)*(N-1) + 1
      DO 10 I = 1,N
         IF (ARRAY(IX) .GT. TARGET) THEN
            NVAL = NVAL + 1
            INDEX(NVAL) = I
         END IF
         IX = IX + INC
   10 CONTINUE
      RETURN
      END
