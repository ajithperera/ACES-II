subroutine write_aces_r2array(mat,dim0,irow,icol,fname,fileout)
   implicit none
   integer dim0,irow,icol,fileout
   character(len=*) fname
   double precision mat(dim0)
   integer i,j,ind,filewrite

   filewrite = 25
   open(unit=filewrite,file=fname,form='formatted',status='unknown')
   ind = 0
   do j = 1, icol
      do i = 1, irow
         ind = ind + 1
         write(filewrite,'(2i5,f20.15)') i,j,mat(ind)
      end do
   end do
   close(filewrite)
end subroutine write_aces_r2array
