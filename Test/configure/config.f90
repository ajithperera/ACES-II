
program config
  integer :: i=0, l, iarr(2)
  integer :: nml0=0, nml1=1
  namelist /struct/ nml0, nml1

  print *, "==> IDACCM test"
  print *, ""
  inquire(iolength=l) i
  print *, "The record length of an integer is ", l
  selectcase(l)
  case(1)
    print *, "If this is a 32-bit build, then RECL is in 4-byte units."
    print *, "Add the define '-D_RECL_IS_WORDS_' to the ACES cpp command."
    print *, "If this is a 64-bit build, then this behavior is unique."
  case(2)
    print *, "If this is a 64-bit build, then RECL is in 4-byte units."
    print *, "Add the define '-D_RECL_IS_WORDS_' to the ACES cpp command."
  case(4)
    print *, "If this is a 32-bit build, then RECL is in 1-byte units."
    print *, "Remove the define '-D_RECL_IS_WORDS_' from the ACES cpp command."
  case(8)
    print *, "If this is a 64-bit build, then RECL is in 1-byte units."
    print *, "Remove the define '-D_RECL_IS_WORDS_' from the ACES cpp command."
  case default
    print *, "Something is wrong with the iolength assignment."
  end select

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
  selectcase(l)
  case(1)
    print *, "The memory address units are in integers."
  case(2)
    print *, "If this is a 64-bit build, &
             &then the memory address units are in 4-byte words."
    print *, "Otherwise, something is wrong."
  case(4)
    print *, "If this is a 32-bit build, &
             &then the memory address units are in bytes."
    print *, "Otherwise, something is wrong."
  case(8)
    print *, "If this is a 64-bit build, &
             &then the memory address units are in bytes."
    print *, "Otherwise, something is wrong."
  case default
    print *, "Something is wrong with the loc function."
  end select

  print *, ""
  print *, "==> namelist formatting test"
  print *, ""
  write(*,nml=struct)

end program config

