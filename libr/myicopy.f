      SUBROUTINE MYICOPY(I1,I2,LEN)
      DIMENSION I1(LEN),I2(LEN)
      DO 10 I=1,LEN
       I2(I)=I1(I)
10    CONTINUE
      RETURN
      END
