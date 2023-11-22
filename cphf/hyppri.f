      SUBROUTINE HYPPRI(AMAT,SPC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION AMAT(3,3,3)
      CHARACTER*(*) SPC
C      FACTOR CONVERTS FROM AU**3 TO ANGSTROM**3
      PARAMETER(FACTOR=0.148 184 D00, ONE=1.0D00)
      PARAMETER(LUCMD=5, LUPRI=6)
      CHARACTER*2 LAB(3),LAB2(6)
C
      IF(INDEX(SPC,'EXP').NE.0) THEN
       FAC=FACTOR
      ELSE
       FAC=ONE
      ENDIF
      IF(INDEX(SPC,'PRI').NE.0) THEN
        LAB(1)='EA'
        LAB(2)='EB'
        LAB(3)='EC'
        WRITE(LUPRI, '(15X,A/)')
     *    '(Along principal axes of moments of inertia)'
      ELSE
        LAB(1)='Ex'
        LAB(2)='Ey'
        LAB(3)='Ez'
        LAB2(1)='xx'
        LAB2(2)='xy'
        LAB2(3)='xz'
        LAB2(4)='yy'
        LAB2(5)='yz'
        LAB2(6)='zz'
      ENDIF
       WRITE(LUPRI,1000) LAB
       WRITE(LUPRI,2000)LAB2(1),FAC*AMAT(1,1,1),FAC*AMAT(1,1,2),
     *                            FAC*AMAT(1,1,3)
       WRITE(LUPRI,2000)LAB2(2),FAC*AMAT(1,1,2),FAC*AMAT(2,1,2),
     *                          FAC*AMAT(3,1,2)
       WRITE(LUPRI,2000)LAB2(3),FAC*AMAT(1,1,3),FAC*AMAT(2,1,3),
     *                           FAC*AMAT(3,1,3)
       WRITE(LUPRI,2000)LAB2(4),FAC*AMAT(1,2,2),FAC*AMAT(2,2,2),
     *                           FAC*AMAT(3,2,2)
       WRITE(LUPRI,2000)LAB2(5),FAC*AMAT(1,2,3),FAC*AMAT(2,2,3),
     *                           FAC*AMAT(3,2,3)
       WRITE(LUPRI,2000)LAB2(6),FAC*AMAT(1,3,3),FAC*AMAT(2,3,3),
     *                           FAC*AMAT(3,3,3)
      WRITE(LUPRI, '(/)')
      RETURN
1000  FORMAT(15X,3(10X,A2)/)
2000  FORMAT(12X,A2,3X,3F12.6)
      END
