      SUBROUTINE DLOOP(NHKT1,NHKT2,KHKT1,KHKT2,DIAG12,IAB0X,IAB0Y,IAB0Z,
     *                 DER0,DER1,DER2,LOOP0,LOOP1,LOOP2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL DIAG12, DER0, DER1, DER2
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN)
      DIMENSION ISTEP(MXAQN), MVALUE(MXAQN), NVALUE(MXAQN)

      DATA ISTEP /1,2,2,3,3,3,4,4,4,4,5,5,5,5,5,
     &            6,6,6,6,6,6,7,7,7,7,7,7,7,
     &            8,8,8,8,8,8,8,8/
      DATA MVALUE /0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,
     &            5,4,3,2,1,0,6,5,4,3,2,1,0,
     &            7,6,5,4,3,2,1,0/
      DATA NVALUE /0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,
     &            0,1,2,3,4,5,0,1,2,3,4,5,6,
     &            0,1,2,3,4,5,6,7/        
C
      IAB1X = IAB0X + 1
      IAB1Y = IAB0Y + 1
      IAB1Z = IAB0Z + 1
      DO 100 ICOMP1 = 1,KHKT1
         LVAL1 = NHKT1 - ISTEP(ICOMP1) + IAB1X
         MVAL1 = MVALUE(ICOMP1) + IAB1Y
         NVAL1 = NVALUE(ICOMP1) + IAB1Z
         IF (DIAG12) THEN
            MAX2 = ICOMP1
         ELSE
            MAX2 = KHKT2
         END IF
         DO 200 ICOMP2 = 1,MAX2
            LVAL12 = LVAL1 + NHKT2 - ISTEP(ICOMP2)
            MVAL12 = MVAL1 + MVALUE(ICOMP2)
            NVAL12 = NVAL1 + NVALUE(ICOMP2)
            NUMT0 = LVAL12/IAB1X
            NUMU0 = MVAL12/IAB1Y
            NUMV0 = NVAL12/IAB1Z
            IF (DER0) LOOP0 = NUMT0*NUMU0*NUMV0
            IF (DER1) THEN
               NUMT1 = (LVAL12 + 1)/IAB1X
               NUMU1 = (MVAL12 + 1)/IAB1Y
               NUMV1 = (NVAL12 + 1)/IAB1Z
               LOOP1 = NUMT1*NUMU0*NUMV0 + NUMT0*NUMU1*NUMV0
     *               + NUMT0*NUMU0*NUMV1
            END IF
            IF (DER2) THEN
               NUMT1 = (LVAL12 + 1)/IAB1X
               NUMU1 = (MVAL12 + 1)/IAB1Y
               NUMV1 = (NVAL12 + 1)/IAB1Z
               NUMT2 = (LVAL12 + 2)/IAB1X
               NUMU2 = (MVAL12 + 2)/IAB1Y
               NUMV2 = (NVAL12 + 2)/IAB1Z
               LOOP2 = NUMT2*NUMU0*NUMV0 + NUMT1*NUMU1*NUMV0
     *               + NUMT1*NUMU0*NUMV1 + NUMT0*NUMU2*NUMV0
     *               + NUMT0*NUMU1*NUMV1 + NUMT0*NUMU0*NUMV2
            END IF
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
