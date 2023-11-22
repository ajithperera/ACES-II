      subroutine prnmat(label,matrix,dim1,dim2,verprn)
      implicit none

      character*20 label
      integer dim1, dim2
      double precision matrix(dim1,dim2)
      logical verprn

      if (verprn) then
         write(*,*)
         write(*,16) label
         call output3(matrix,1,dim1,1,dim2,dim1,dim2,1)
         write(*,*)
      end if

 16   format(A20,' Matrix')

      return
      end
