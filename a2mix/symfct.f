      SUBROUTINE SYMFCT(ISYMYZ,ISYMXZ,ISYMXY,ISYMRX,ISYMRY,ISYMRZ,ISYMI,
     &   NATOMS,NCNTR,COORD,IQ1,IQ2,IQ3,IQ4,IQ5,IQ6,IQ7,IQ8,IFACT)
C
C     This routine determines which quadrants to integrate over and the 
C     factor to multiply the symmetry integration by.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL C1LT,C2LT,C3LT
C
      PARAMETER (EPS = 1.0D-10)
C
      DIMENSION COORD(3*NATOMS)
C
      COMMON /IPAR/ LUOUT
C
C     Quadrants are labeled
C     xyz=1, xy(-z)=2, x(-y)z=3, x(-y)(-z)=4
C     (-x)yx=5, (-x)y(-z)=6, (-x)(-y)z=7, (-x)(-y)(-z)=8
C 
C     Determine which symmetry elements apply to atom NCNTR
C
      IYZ=ISYMYZ
      IXZ=ISYMXZ
      IXY=ISYMXY
      IRX=ISYMRX
      IRY=ISYMRY
      IRZ=ISYMRZ
      IMI=ISYMI
C
C     Determine which coordinates are less than EPS away from zero.
C
      C1LT = ABS(COORD(3*(NCNTR-1)+1)) .LT. EPS
      C2LT = ABS(COORD(3*(NCNTR-1)+2)) .LT. EPS
      C3LT = ABS(COORD(3*(NCNTR-1)+3)) .LT. EPS
C
C     Set zeros for COORD matrix
c      DO 999 IV=1,3*NATOMS
c         write(*,*) 'coord',iv,'=',coord(iv)
c         IF(DABS(COORD(IV)).LT.1.D-10) COORD(IV)=0.D+00
c  999 CONTINUE
      IF (.NOT.C1LT .AND. .NOT.C2LT .AND. .NOT.C3LT) THEN
        ISYMYZ=0
        ISYMXZ=0
        ISYMXY=0
        ISYMRX=0
        ISYMRY=0
        ISYMRZ=0
        ISYMI=0
C
      ELSE IF (C1LT .AND. .NOT.C2LT .AND. .NOT.C3LT) THEN
        ISYMXZ=0
        ISYMXY=0
        ISYMRX=0
        ISYMRY=0
        ISYMRZ=0
        ISYMI=0
C
      ELSE IF (.NOT.C1LT .AND. C2LT .AND. .NOT.C3LT) THEN
        ISYMYZ=0
        ISYMXY=0
        ISYMRX=0
        ISYMRY=0
        ISYMRZ=0
        ISYMI=0
C
      ELSE IF (.NOT.C1LT .AND. .NOT.C2LT .AND. C3LT) THEN
        ISYMYZ=0
        ISYMXZ=0
        ISYMRX=0
        ISYMRY=0
        ISYMRZ=0
        ISYMI=0
C
      ELSE IF (C1LT .AND. C2LT .AND. .NOT.C3LT) THEN
        ISYMXY=0
        ISYMRX=0
        ISYMRY=0
        ISYMI=0
C
      ELSE IF (C1LT .AND. .NOT.C2LT .AND. C3LT) THEN
        ISYMXZ=0
        ISYMRX=0
        ISYMRZ=0
        ISYMI=0
C
      ELSE IF (.NOT.C1LT .AND. C2LT .AND. C3LT) THEN
        ISYMYZ=0
        ISYMRY=0
        ISYMRZ=0
        ISYMI=0
      ENDIF
C
C     Initial symmetry factor
      IFACT=1
C
C     Determine which quadrants are equivalent
      IQ1=1
      IQ2=1
      IQ3=1
      IQ4=1
      IQ5=1
      IQ6=1
      IQ7=1
      IQ8=1
C
      IF (ISYMYZ.EQ.1.OR.ISYMRY.EQ.1.OR.ISYMRZ.EQ.1.OR.ISYMI.EQ.1) THEN
        IQ5=0
        IQ6=0
        IQ7=0
        IQ8=0
        IF (ISYMXZ.EQ.1.OR.ISYMRX.EQ.1) THEN
          IQ3=0
          IQ4=0
          IF (ISYMXY.EQ.1) THEN
            IQ2=0
            IFACT=8
          ELSE
            IFACT=4
          ENDIF
        ELSE IF (ISYMXY.EQ.1) THEN
          IQ2=0
          IQ4=0
          IFACT=4
        ELSE
          IFACT=2
        ENDIF
      ELSE IF (ISYMXZ.EQ.1.OR.ISYMRX.EQ.1) THEN
        IQ3=0
        IQ4=0
        IQ7=0
        IQ8=0
        IF (ISYMXY.EQ.1) THEN
          IQ2=0
          IQ6=0
          IFACT=4
        ELSE
          IFACT=2
        ENDIF
      ELSE IF (ISYMXY.EQ.1) THEN
        IQ2=0
        IQ4=0
        IQ6=0
        IQ8=0
        IFACT=2
      ENDIF
C
 300  CONTINUE
C
      ISYMYZ=IYZ
      ISYMXZ=IXZ
      ISYMXY=IXY
      ISYMRX=IRX
      ISYMRY=IRY
      ISYMRZ=IRZ
      ISYMI=IMI
C
      RETURN
      END
