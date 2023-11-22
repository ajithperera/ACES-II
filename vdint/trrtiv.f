      SUBROUTINE TRRTIV(GM,X,Y,Z,NUCDEP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C  THIS SUBROUTINE CHECKS IF THE MOLECULAR GRADIENT
C  SATISFY THE TRANLATIONAL-ROTATIONAL INVARIANCE 
C  CONDITIONS
C
C  J. GAUSS , NEVER TO BE PUBLISHED (IN CONTRAST TO
C                                    TUH'S FORMULAS)
C
C  THE OTHER DIFFERENCE IS THAT THIS ROUTINE
C  WORKS WHILE TUH'S ORIGINAL ROUTINE FAILED
C  QUITE OFTEN BADLY.   WELL, AND THIS GOT PUBLISHED !!!
C
CEND
C
C  WRITTEN JAN/91  JG
C
      DIMENSION GM(3,1),X(1),Y(1),Z(1)
      COMMON/IPRI/IPRINT
      DATA ZERO,TRESH1,TRESH2/0.D0,1.D-08,1.D-05/
C
      IF(IPRINT.NE.0) THEN
       WRITE(6,900)
900    FORMAT('  @TRRTIV-I, Check for translational-rotational',
     &       ' invariances.')
       WRITE(6,901)
901    FORMAT(' ',72('-'))
      ENDIF
      REWIND(UNIT=81)
C
C  READ IN CARTESIAN COORDINATES OF ATOMS
C
      READ(81,3000) NUCDEP,ENERGY
3000  FORMAT(I5,F20.10)
C
      DO 10 IATOM=1,NUCDEP
       READ(81,4000) XJUNK,X(IATOM),Y(IATOM),Z(IATOM)
10    CONTINUE
      DO 20 IATOM=1,NUCDEP
       READ(81,4000) XJUNK,GM(1,IATOM),GM(2,IATOM),GM(3,IATOM)
20    CONTINUE
C
4000  FORMAT(4F20.10)
C
C  CHECK FIRST FOR TRANLATIONAL INVARIANCE
C
      TRX=ZERO
      TRY=ZERO
      TRZ=ZERO
      DO 100 IATOM=1,NUCDEP
       TRX=TRX+GM(1,IATOM)
       TRY=TRY+GM(2,IATOM)
       TRZ=TRZ+GM(3,IATOM)
100   CONTINUE 
C
C  WRITE OUT RESULT FOR TRANSLATIONAL INVARIANCE
C
      IF(MAX(ABS(TRX),ABS(TRY),ABS(TRZ)).LE.TRESH1) THEN
       IF(IPRINT.NE.0) THEN
        write(6,5000)
5000    FORMAT('  Molecular gradient satisfies translational',
     &         ' invariance condition.')
        write(6,5001) TRX
        write(6,5002) TRY
        write(6,5003) TRZ
5001    FORMAT('   translational invariance in x direction :   ',E10.4)
5002    FORMAT('   translational invariance in y direction :   ',E10.4)
5003    FORMAT('   translational invariance in z direction :   ',E10.4)
       ENDIF
      ELSE
       write(6,5010)
5010   FORMAT('  Molecular gradient does not satisfy translational',
     &        ' invariance condition.')
       write(6,5001) TRX
       write(6,5002) TRY
       write(6,5003) TRZ
C      CALL ERREX
      ENDIF
C
C  CHECK ROTATIONAL INVARIANCE
C
      RRX=ZERO
      RRY=ZERO
      RRZ=ZERO
      DO 200 IATOM=1,NUCDEP
       RRX=RRX+GM(2,IATOM)*Z(IATOM)-GM(3,IATOM)*Y(IATOM)
       RRY=RRY+GM(3,IATOM)*X(IATOM)-GM(1,IATOM)*Z(IATOM)
       RRZ=RRZ+GM(1,IATOM)*Y(IATOM)-GM(2,IATOM)*X(IATOM)
200   CONTINUE
C
C  WRITE OUT RESULT FOR ROTATIONAL INVARIANCE
C
      IF(MAX(ABS(RRX),ABS(RRY),ABS(RRZ)).LE.TRESH2) THEN
       IF(IPRINT.NE.0) THEN
        write(6,6000)
6000    FORMAT('  Molecular gradient satisfies rotational',
     &         ' invariance condition.')
        write(6,6001) RRX
        write(6,6002) RRY
        write(6,6003) RRZ
6001    FORMAT('   rotational invariance in yz-plane :        ',E10.4)
6002    FORMAT('   rotational invariance in xz-plane :        ',E10.4)
6003    FORMAT('   rotational invariance in xy-plane :        ',E10.4)
       ENDIF
      ELSE
       write(6,6010)
6010   FORMAT('  Molecular gradient does not satisfy rotational',
     &        ' invariance condition.')
       write(6,6001) RRX
       write(6,6002) RRY
       write(6,6003) RRZ
C      CALL ERREX
      ENDIF
C
      IF(IPRINT.NE.0) write(6,901)
C
      CLOSE(UNIT=81,STATUS='KEEP')
C
      RETURN
      END
