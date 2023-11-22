
      program config
      implicit none

      integer l, iarr(2)
      integer nml0, nml1
      namelist /struct/ nml0, nml1

      print *, "==> IDACCM test"
      print *, ""
      print *, "This test requires a Fortran 90 compiler."

      print *, ""
      print *, "==> character translation intrinsics test"
      print *, ""
      print *, "iachar('A') = ", iachar('A'), " (should be 65)"
      print *, " ichar('A') = ",  ichar('A'), " ( could be 65)"
      print *, " achar( 65) = '", achar(65), "' (should be 'A')"
      print *, "  char( 65) = '",  char(65), "' ( could be 'A')"

      print *, ""
      print *, "==> memory address resolution test"
      print *, ""
      l=loc(iarr(2))-loc(iarr(1))
      print *, "Two consecutive integers are ",l," memory unit(s) apart"
      if      (l.eq.1) then
         print *, "The memory address units are in integers."
      else if (l.eq.2) then
         print *, "If this is a 64-bit build, ",
     &            "then the memory address units are in 4-byte words."
         print *, "Otherwise, something is wrong."
      else if (l.eq.4) then
         print *, "If this is a 32-bit build, ",
     &            "then the memory address units are in bytes."
         print *, "Otherwise, something is wrong."
      else if (l.eq.8) then
         print *, "If this is a 64-bit build, ",
     &            "then the memory address units are in bytes."
         print *, "Otherwise, something is wrong."
      else
         print *, "Something is wrong with the loc function."
      end if

      print *, ""
      print *, "==> namelist formatting test"
      print *, ""
      nml0=0
      nml1=1
      write(*,nml=struct)

      end

