
C Subroutine to evaluate fourth-order energy from tripleS
c amplitudes and denominators for given ijk (l Signifies 
C we are computing 4th order contribution to L-CCSD(T)

c#define _DEBUG_E4TIJKL

      SUBROUTINE E4TIJKL(T3,L3,D3,LENGTH,DIJK,E4TRIP,caller)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION T3(LENGTH), D3(LENGTH), L3(LENGTH)
      character*(*) caller
      INTEGER ABC
#ifdef _DEBUG_E4TIJKL
      integer iLength, len
      iLength = len(caller)
#endif /* _DEBUG_E4TIJKL */
      DO ABC = 1, LENGTH
         E4TRIP = E4TRIP + (DIJK - D3(ABC)) * T3(ABC) * L3(ABC)
#ifdef _DEBUG_E4TIJKL
c         if (caller(1:iLength).eq.'trps3.F') then
            write(*,*) caller(1:iLength),' : e4trip += ',
     &                 (dijk - d3(abc)),' * (',t3(abc),')^2'
c         end if
#endif /* _DEBUG_E4TIJKL */
      END DO
      RETURN
      END

