
      subroutine shifts

      integer ishft, rshift, lshift
      intrinsic ishft, rshift, lshift

      print *, "==> bit shifting intrinsic functions test"
      print *, ""
      print *, "rshift(2, 1) = ", rshift(2, 1), " (should be 1)"
      print *, " ishft(2,-1) = ",  ishft(2,-1), " (should be 1)"
      print *, " ishft(2, 1) = ",  ishft(2, 1), " (should be 4)"
      print *, "lshift(2, 1) = ", lshift(2, 1), " (should be 4)"
      print *, ""

      return
      end

