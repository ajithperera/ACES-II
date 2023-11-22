
c mangling table: (have to use doubles since int size changes)
c  1.d0 = normal	(sub -> sub)
c  2.d0 = underscore	(sub -> sub_)
c  3.d0 = uppercase	(sub -> SUB)

      subroutine fortran_names
      double precision d
      common /ftncom/ d

      print *, "==> FORTRAN name-mangling test"
      print *, ""
      d=0.d0
      call c_sub
      if (0.5d0.lt.d.and.d.lt.1.5d0) then
         print *, "Common block names are direct translation."
         print *, "> Remove '-DCB_SUFFIX -DCB_UPPER'"
      else
         if (1.5d0.lt.d.and.d.lt.2.5d0) then
            print *, "Common block names are suffixed."
            print *, "> Add '-DCB_SUFFIX', remove '-DCB_UPPER'"
         else
            if (2.5d0.lt.d.and.d.lt.3.5d0) then
               print *, "Common block names are uppercase."
               print *, "> Add '-DCB_UPPER', remove '-DCB_SUFFIX'"
            else
               print *, "Something is wrong."
               print *, "> ERROR: d=",d
            end if
         end if
      end if
      print *, ""

      return
      end

