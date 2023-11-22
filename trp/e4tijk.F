
C SUBROUTINE TO EVALUATE FOURTH-ORDER ENERGY FROM TRIPLES
C AMPLITUDES AND DENOMINATORS FOR GIVEN IJK.

c#define _DEBUG_E4TIJK

      SUBROUTINE E4TIJK(T3,D3,LENGTH,DIJK,E4TRIP,caller)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION T3(LENGTH), D3(LENGTH)
      character*(*) caller
      INTEGER ABC
#ifdef _DEBUG_E4TIJK
      integer iLength, len
      iLength = len(caller)
#endif /* _DEBUG_E4TIJK */
      DO ABC = 1, LENGTH
         E4TRIP = E4TRIP + (DIJK - D3(ABC)) * T3(ABC) * T3(ABC)
#ifdef _DEBUG_E4TIJK
c         if (caller(1:iLength).eq.'trps3.F') then
            write(*,*) caller(1:iLength),' : e4trip += ',
     &                 (dijk - d3(abc)),' * (',t3(abc),')^2'
c         end if
#endif /* _DEBUG_E4TIJK */
      END DO
      RETURN
      END

