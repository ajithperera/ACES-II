
c This function returns the address of the first element in array(n*inc)
c that is not equal to target.
c BEWARE - ISRCHNE is not defined if n=0

      INTEGER FUNCTION ISRCHNE(N,ARRAY,INC,TARGET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER ARRAY,TARGET
      DIMENSION ARRAY(N)
#ifdef _STRINGENT_CHECK
      if (n.lt.0)   dtmp = dsqrt(dfloat(n))
      if (n.eq.0)   dtmp = 1/n
      if (inc.lt.0) dtmp = dsqrt(dfloat(inc))
      if (inc.eq.0) dtmp = 1/inc
#endif /* _STRINGENT_CHECK */
      ISRCHNE = 1
      IF (INC.LT.0) THEN
         JLOOK = 1 + INC + ( -N * INC )
         INCR = -INC
      ELSE
         JLOOK = 1
         INCR = INC
      END IF
      DO I = 1, N
         IF (ARRAY(JLOOK).NE.TARGET) GOTO 200
         JLOOK = JLOOK + INC
         ISRCHNE = ISRCHNE + INCR
      END DO
200   CONTINUE
      RETURN
      END
