      SUBROUTINE ADDONEH(W,SCR,ISIZE,FVV,FOO,DISSZ1,NUMDS1,ISPIN)
C
C THIS SUBROUTINE ACCEPTS THE QUANTITIES F(AB) [FVV], F(IJ) [FOO]
C AND W(AB,IJ) AND RETURNS
C
C   Z(AB,IJ) = W(AB,IJ) + F(AB) DELTA(IJ) - F(IJ) DELTA (AB)
C
C WHERE DELTA IS THE KRONECKER DELTA FUNCTION.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSZ1,POP,VRT
      DIMENSION W(ISIZE),SCR(ISIZE),FOO(*),FVV(*)
C
      COMMON/SYMINF/NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
C 
      DATA ONE,ONEM /1.0D0,-1.0D0/
C
      IOFFH0=1
      DO 20 IRREP=1,NIRREP
       NUMO=POP(IRREP,ISPIN)
       IOFFHO=IOFFH0
       DO 21 INDXJ=1,POP(IRREP,ISPIN)
        CALL SAXPY(NFEA(ISPIN),ONE,FVV,1,W(IOFFHO),1)
        IOFFHO=IOFFHO+DISSZ1*(NUMO+1)
21     CONTINUE
       IOFFH0=IOFFH0+NUMO*NUMO*DISSZ1
20    CONTINUE
      CALL TRANSP(W,SCR,NUMDS1,DISSZ1)
      IOFFH0=1
      DO 22 IRREP=1,NIRREP
       NUMV=VRT(IRREP,ISPIN)
       IOFFHO=IOFFH0
       DO 23 INDXB=1,VRT(IRREP,ISPIN)
        CALL SAXPY(NFMI(ISPIN),ONEM,FOO,1,SCR(IOFFHO),1)
        IOFFHO=IOFFHO+NUMDS1*(NUMV+1)
23     CONTINUE
       IOFFH0=IOFFH0+NUMV*NUMV*NUMDS1
22    CONTINUE
      CALL TRANSP(SCR,W,DISSZ1,NUMDS1)
      RETURN
      END