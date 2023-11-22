      SUBROUTINE LMNVAL(NHKTA,IBCMPA,LVALUE,MVALUE,NVALUE)
C
C     This subroutine determines LVALUE(KHKTA), MVALUE(KHKTA), and
C     NVALUE(KHKTA) from NHKTA and IBCMPA.
C
C     tuh March 87
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2)
      DIMENSION LVALUE(MXAQN), MVALUE(MXAQN), NVALUE(MXAQN),
     *          ISTEP(MXAQN),  MVAL(MXAQN),   NVAL(MXAQN)

      DATA ISTEP /1,2,2,3,3,3,4,4,4,4,5,5,5,5,5,
     *            6,6,6,6,6,6,7,7,7,7,7,7,7,
     *            8,8,8,8,8,8,8,8/ 

      DATA MVAL  /0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2,1,0,
     *            6,5,4,3,2,1,0,7,6,5,4,3,2,1,0/

      DATA NVAL  /0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,1,2,3,4,5,
     *            0,1,2,3,4,5,6,0,1,2,3,4,5,6,7/

      IBTSHR(I,J) = ISHFT(I,-J)
      IBTAND(I,J) = AND(I,J)
C
      ICOMP = 0
      DO 100 I = 1, NHKTA*(NHKTA+1)/2
         IF (IBTAND(IBTSHR(IBCMPA,I-1),1) .EQ. 1) THEN
            ICOMP = ICOMP + 1
            LVALUE(ICOMP) = NHKTA - ISTEP(I)
            MVALUE(ICOMP) = MVAL(I)
            NVALUE(ICOMP) = NVAL(I)
         END IF
  100 CONTINUE
      RETURN
      END
