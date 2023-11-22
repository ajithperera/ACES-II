      SUBROUTINE SXF2(IRREPX,FDAI,FIJ,SDAI,POP,VRT)
C
C THIS ROUTINE CALCULATES THE FOLLOWING CONTRIBUTIONS
C TO THE ROHF-FOCK MATRIX DERIVATIVES
C
C   - SUM J S(A,J)^chi F(J,I)
C
CEND
C
C CODED APRIL/91 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
C
      DIMENSION FDAI(1),FIJ(1),SDAI(1),POP(8),VRT(8)
C
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
      DATA ONE,ONEM,HALFM /1.D0,-1.0D0,-0.5D0/ 
C
C LOOP OVER ALL IRREPS
C
      IOFFFD=1
      IOFFR=1
C
      DO 1000 IRREPR=1,NIRREP
C
       IRREPL=DIRPRD(IRREPX,IRREPR)
C
       NUMR=POP(IRREPR)
       NUML=VRT(IRREPL)
C
C MULTIPLY BY RIGHT
C
       CALL XGEMM('N','N',NUML,NUMR,NUMR,ONEM,
     &            SDAI(IOFFFD),NUML,FIJ(IOFFR),NUMR,
     &            ONE,FDAI(IOFFFD),NUML)
C
       IOFFFD=IOFFFD+NUML*NUMR
       IOFFR=IOFFR+NUMR*NUMR
C
1000  CONTINUE
C
      RETURN
      END
