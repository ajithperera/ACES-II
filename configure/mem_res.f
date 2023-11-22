
c int iIntLn = byte-length of the default integer type

      subroutine mem_res(iIntLn)
      integer iIntLn, l, iArr(2)

      print *, "==> memory address resolution test"
      print *, ""
      l=loc(iArr(2))-loc(iArr(1))
      print *, "Two consecutive integers are ",l," unit(s) apart."
      if (l.eq.iIntLn) then
         print *, "LOC is in bytes."
         print *, "> Remove '-D_PTRS_ARE_WORDS'"
      else
         if (l.eq.1) then
            print *, "LOC is in integer words."
            print *, "> Add '-D_PTRS_ARE_WORDS'"
         else
            print *, "Something is wrong."
            print *, "> ERROR: iIntLn=",iIntLn,", loc2-loc1=",l
         end if
      end if
      print *, ""

      return
      end

