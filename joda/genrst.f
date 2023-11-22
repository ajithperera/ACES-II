       SUBROUTINE GENRST(J,SCR,NTIME,NUMTIM,IORDGP,JX,CLSTYP)
C
C THIS SUBROUTINE GENERATES THE REMAINING TRANSFORMATION MATRICES
C  BELONGING TO A GIVEN CLASS BY SUCCESSIVELY APPLYING A PARTICULAR
C  TRANSFORMATION MATRIX.  PARTICULARLY USEFUL FOR NON-CUBIC GROUPS.
C
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION SCR(9*IORDGP)
       INTEGER CLSTYP(IORDGP) 
       IP(I)=1+9*(I-1) 
       DO 10 I=1,NUMTIM
        NTIME=NTIME+1
        CALL UNITRY(SCR(IP(J)),SCR(IP(NTIME-1)),SCR(IP(NTIME)),3)
        CLSTYP(NTIME)=JX
10     CONTINUE
       RETURN
       END