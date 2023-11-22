
! int iIntLn = byte-length of the default integer type

subroutine recl_units(iIntLn)
   integer :: iIntLn, l

   print *, "==> RECL units test"
   print *, ""
   inquire(iolength=l) iIntLn
   print *, "The record length of an integer is ",l," units."
   if (iIntLn==l) then
      print *, "RECL is in bytes."
      print *, "> Remove '-D_RECL_IS_WORDS_'"
   else
      if (iIntLn==4*l) then
         print *, "RECL is in 4-byte words."
         print *, "> Add '-D_RECL_IS_WORDS_'"
      else
         print *, "Something is wrong."
         print *, "> ERROR: iIntLn=",iIntLn,", iolength=",l
      end if
   end if
   print *, ""

   return
end

