      SUBROUTINE MSZPRI(AMAT2,AMAT1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL PARA
      DIMENSION AMAT2(3,3),AMAT1(3,3)
C      FACTOR CONVERTS FROM AU**3 TO ANGSTROM**3
      PARAMETER(FACTOR=0.148 184 D0)
      CHARACTER*2 LAB(3)
C
      DATA ONE,THREE,FOUR,FOURTH /1.D0,3.D0,4.D0,0.25D0/
C
      CALL HEADER(
     & 'Paramagnetic MBPT(2) contribution of magnetizability tensor',
     &  -1)
C
      LAB(1)='Bx'
      LAB(2)='By'
      LAB(3)='Bz'
      WRITE(*, 1000) LAB
      WRITE(*, 2000) LAB(1),AMAT2(1,1),AMAT2(1,2),
     &                           AMAT2(1,3)
      WRITE(*, 2000) LAB(2),AMAT2(2,1),AMAT2(2,2),
     &                           AMAT2(2,3)
      WRITE(*, 2000) LAB(3),AMAT2(3,1),AMAT2(3,2),
     &                           AMAT2(3,3)
      WRITE(*, '(/)')
C
C SAVE MAGNETIZABILITY TENSOR ON FILE CHI
C
      OPEN(UNIT=82,FILE='CHI',STATUS='UNKNOWN',FORM='FORMATTED')
      REWIND(82)
      READ(82,3000) (AMAT1(I,1),I=1,3)
      READ(82,3000) (AMAT1(I,2),I=1,3)
      READ(82,3000) (AMAT1(I,3),I=1,3)
C
      CALL SAXPY(9,ONE,AMAT1,1,AMAT2,1)
C
      REWIND(82)
      WRITE(82,3000) (AMAT2(I,1),I=1,3)     
      WRITE(82,3000) (AMAT2(I,2),I=1,3)     
      WRITE(82,3000) (AMAT2(I,3),I=1,3)
C
      CLOSE(UNIT=82,STATUS='KEEP')
C
      CALL HEADER('Total MBPT(2) magnetizability tensor',-1)
      WRITE(*,1000) LAB
      WRITE(*,2000) LAB(1),AMAT2(1,1),AMAT2(1,2),
     &                         AMAT2(1,3)
      WRITE(*,2000) LAB(2),AMAT2(2,1),AMAT2(2,2),
     &                         AMAT2(2,3)
      WRITE(*,2000) LAB(3),AMAT2(3,1),AMAT2(3,2),
     &                         AMAT2(3,3)
      WRITE(*, '(/)')
C
      AMAG=(AMAT2(1,1)+AMAT2(2,2)+AMAT2(3,3))/THREE
C
      WRITE(*,4000) AMAG
      WRITE(*, '(/)')
C
      RETURN
1000  FORMAT(15X,3(10X,A2)/)
2000  FORMAT(12X,A2,3X,3F12.6)
3000  FORMAT(3F20.10)
4000  FORMAT(T15,'Total MBPT(2) magnetizability :',F10.4)
      END
