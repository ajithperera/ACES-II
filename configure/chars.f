
      subroutine chars

      print *, "==> character intrinsic functions test"
      print *, ""
      print *, "iachar('A') = ", iachar('A'), " (should be 65)"
      print *, " ichar('A') = ",  ichar('A'), " ( could be 65)"
      print *, " achar( 65) = '", achar(65), "' (should be 'A')"
      print *, "  char( 65) = '",  char(65), "' ( could be 'A')"
      print *, ""

      return
      end

