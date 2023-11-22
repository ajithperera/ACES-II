C
C ***************************************************************
C  SUBROUTINES THAT RELATE TO NORMALIZATION OF S-VECTORS
C ***************************************************************
C
      SUBROUTINE NORMVEC(VEC, NSIZEC, SCR, MAXCOR, IUHF, ICALC, IRREP)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION VEC(NSIZEC), SCR(MAXCOR)
C
      I010 = 1
      I020 = I010 + NSIZEC
      CALL SCOPY(NSIZEC,VEC,1,SCR(I010),1)
      IF(IUHF.EQ.0)THEN
         CALL SPNTSING(NSIZEC,SCR(I010),SCR(I020),MAXCOR-I020+1,
     $      IRREP,ICALC)
      ENDIF
      Z=SDOT(NSIZEC,VEC,1,SCR(I010),1)
      X=1.0D0/SQRT(Z)
      CALL SSCAL (NSIZEC,X,VEC,1)
C
      RETURN
      END
