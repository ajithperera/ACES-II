      SUBROUTINE LARM(I1,I2)
      INTEGER I1(50),I2(20)
      WRITE (6,*) '@LARM: NOT ENOUGH MEMORY!!!!!'
      WRITE (6,*) 'I1(1:50) = ',I1
      WRITE (6,*) 'I2(1:20) = ',I2
      CALL ERREX
      END
