      SUBROUTINE CLOSMO_C4(IUNIT)
      COMMON/FSTRIPE/NUMSPLIT(5)
      nsplit=numsplit(iunit-49)
      do 333 i=1,nsplit
      CLOSE(UNIT=IUNIT+5*(i-1),STATUS='KEEP')
333   continue
      RETURN
      END
