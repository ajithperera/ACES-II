
      subroutine args
      implicit none
      integer i
      integer*4 i4
      character*72 arg(2)

      print *, "==> getarg integer argument test"
      print *, ""
      do i = 0, 1
         call getarg(i,arg(1+i))
      end do
      if (arg(1).eq.arg(2)) then
         print *, "getarg does not use the default integer."
         do i4 = 0, 1
            call getarg(i4,arg(1+i4))
         end do
         if (arg(1).eq.arg(2)) then
            print *, "Something is wrong."
            print *, "> ERROR: args are"
            print *, "> 0: '",arg(1),"'"
            print *, "> 1: '",arg(2),"'"
         else
            print *, "getarg uses 32-bit integers."
            print *, "> Define INTTYPE in tools/f_getarg as 'integer*4'"
         end if
      else
         print *, "getarg uses the default integer."
         print *, "> Define INTTYPE in tools/f_getarg as 'integer'"
      end if
      print *, ""

      return
      end

