      SUBROUTINE E4TIJK(T3,D3,LEN,DIJK,E4TRIP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER ABC
      DIMENSION T3(LEN),D3(LEN)
C
C     SUBROUTINE TO EVALUATE FOURTH-ORDER ENERGY FROM TRIPLES
C     AMPLITUDES AND DENOMINATORS FOR GIVEN IJK.
C
      DO   10 ABC = 1,LEN
      E4TRIP = E4TRIP + (DIJK - D3(ABC)) * T3(ABC) * T3(ABC)
   10 CONTINUE
      RETURN
      END
