C
      SUBROUTINE DRVDSO (DENS, DSOINT, SAVE, NMO, NAO, NCENTR, IUHF)
C 
C Compute the orbital diamagnetic contribution to the NMR  
C spin-spin coupling constant. Coded Ajith 08/93          
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION DENS(NAO*NAO), DSOINT(NAO*NAO), DSO(3),
     &          SAVE(3*NCENTR, 3*NCENTR)
      CHARACTER*2 LABELS(3)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
C
      DATA LABELS /'XX','YY','ZZ'/
C
      NNM1O2(IX) = (IX*(IX-1))/2
      IEXTI(IX)  = 1 + (-1+INT(DSQRT(8.D0*IX+0.999D0)))/2
      IEXTJ(IX)  = IX - NNM1O2(IEXTI(IX))
C
      IMXATM = 50
      IONE   = 1
      INDEX  = 1
      IATOM  = 0
      IRWND  = 0
C
      CALL ZERO(DSOINT, NAO*NAO)
      CALL ZERO(DSO, 3)
      CALL ZERO(SAVE, 9*NCENTR*NCENTR)
C 
      DO 10 ICENTR = 1, NCENTR
         INDEX = INDEX + 1
         DO 20 JCENTR = INDEX, NCENTR
            CALL ZERO(DSO, 3)
            CALL SEEKLB('  ODXX  ', IERR, IRWND, 30)
            IF(IERR .NE. 0) RETURN
            IRWND = 1
            IATOM = IATOM + 1
C
            CALL COMPPR (DSO(1), DENS, DSOINT, NAO, IUHF)
            CALL SEEKLB ('  ODYY  ', IERR, IRWND, 30)
            CALL COMPPR (DSO(2), DENS, DSOINT, NAO, IUHF)
            CALL SEEKLB ('  ODZZ  ',IERR,IRWND, 30)
            CALL COMPPR (DSO(3), DENS, DSOINT, NAO, IUHF)
C
            DO 30 KK = 1, 3
               II = (ICENTR - 1)*3 + KK
               JJ = (JCENTR - 1)*3 + KK
               SAVE(II, JJ) = DSO(KK)
               SAVE(JJ, II) = DSO(KK)
 30         CONTINUE
C
 20      CONTINUE
 10   CONTINUE
C
      RETURN
      END 
