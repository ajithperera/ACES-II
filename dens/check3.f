      SUBROUTINE CHECK3
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/ADD/SUM
      return
      open(unit=76,file='KNUJ',status='UNKNOWN',form='FORMATTED')
      write(76,*) sum
      return
      end
