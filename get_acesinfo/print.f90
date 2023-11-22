subroutine print_rmat(mat,dim0,irow,icol,fileout,itype)
   implicit none
   integer dim0,irow,icol,fileout,itype
   integer i,j,k,ind,rest,num2,iline,istart,iend
   double precision mat(dim0),thres
   double precision, allocatable :: mat2(:,:)

   iline = 5
   thres = 1.0D-7

   if (itype==1) then
      do i = 1, irow
         write (fileout,'(i7,f15.7)') i, mat(i)
      end do
   else if (itype.ge.2) then
      allocate (mat2(irow,icol))
      ind = 0
      do j = 1, icol
         do i = 1, irow
            ind = ind + 1
            mat2(i,j) = mat(ind)
         end do
      end do

      ! icol = iline*num2 + rest
      rest = mod(icol,iline)
      num2 = (icol-rest)/iline
      do k = 1, num2
         istart = iline*(k-1) + 1
         iend = iline*(k-1) + iline
         write (fileout,'(/,7x,10(5x,i5,5x))') (ind,ind=istart,iend)
         if (itype==2) then
            !print full rows
            do i = 1, irow
               write (fileout,'(i7,10f15.7)') i,(mat2(i,ind),ind=istart,iend)
            end do
         elseif (itype==3) then
            !print non-zero rows
            do i = 1, irow
               do j = istart, iend
                  if (abs(mat2(i,j))>=thres) then
                     write (fileout,'(i5,10f15.7)') i,(mat2(i,ind),ind=istart,iend)
                     exit
                  endif
               enddo
            enddo
         endif
      enddo

      if (rest>0) then
         istart = iline*num2+1
         iend = iline*num2+rest
         write (fileout,'(/,7x,10(5x,i5,5x))') (ind,ind=istart,iend)
         if (itype==2) then
            !print full rows
            do i = 1, irow
               write (fileout,'(i5,10f15.7)') i,(mat2(i,ind),ind=istart,iend)
            enddo
         elseif (itype==3) then
            !print non-zero rows
            do i = 1, irow
               do j = istart, iend
                  if (abs(mat2(i,j))>=thres) then
                     write (fileout,'(i5,10f15.7)') i,(mat2(i,ind),ind=istart,iend)
                     exit
                  endif
               enddo
            enddo
         endif
      endif
      !deallocate (mat2)
   endif
end subroutine print_rmat


subroutine print_r2mat_tri(mat1,dim0,dim1,iline,fileout,itype)
   implicit none
   integer dim0,dim1,iline,fileout,itype,ind,i,j
   double precision mat1(dim1) 
   double precision mat2(dim0,dim0) 

   write(fileout,'(a)') 'array'
   do i = 1,dim1
      write(fileout,*) i,mat1(i)
   enddo

   write(fileout,'(a)') 'mat'

   mat2(1:dim0,1:dim0)=0.0d0
   ind=1
   do i = 1,dim0
      do j = 1,i
!        mat2(i,j)=mat1(ind)
         mat2(j,i)=mat1(ind)
         ind=ind+1
      enddo
   enddo
   call print_r2mat(mat2,dim0,dim0,5,fileout,1)
end subroutine print_r2mat_tri


subroutine print_r2mat(mat2,irow,icol,iline,fileout,itype)
   implicit none
   integer irow,icol,fileout,itype
   integer i,j,k,ind,rest,num2,iline,istart,iend
   double precision mat2(irow,icol),thres

   thres = 1.0D-7
   ! icol = iline*num2 + rest
   rest = mod(icol,iline)
   num2 = (icol-rest)/iline
   do k = 1, num2
      istart = iline*(k-1) + 1
      iend = iline*(k-1) + iline
      write (fileout,'(/,5x,10(5x,i5,5x))') (ind,ind=istart,iend)

      if (itype==1) then
         !print full rows
         do i = 1, irow
            write (fileout,'(i5,10f15.7)') i,(mat2(i,ind),ind=istart,iend)
         end do
      elseif (itype==2) then
         !print non-zero rows
         do i = 1, irow
            do j = istart, iend
               if (abs(mat2(i,j))>=thres) then
                  write (fileout,'(i5,10f15.7)') i,(mat2(i,ind),ind=istart,iend)
                  exit
               endif
            enddo
         end do
      else
         write (fileout,*) 'Error: unknown itype=',itype
         exit
      endif
   end do

   if (rest>0) then
      istart = iline*num2+1
      iend = iline*num2+rest
      write (fileout,'(/,5x,10(5x,i5,5x))') (ind,ind=istart,iend)
      if (itype==1) then
         !print full rows
         do i = 1, irow
            write (fileout,'(i5,10f15.7)') i,(mat2(i,ind),ind=istart,iend)
         enddo
      elseif (itype==2) then
         !print non-zero rows
         do i = 1, irow
            do j = istart, iend
               if (abs(mat2(i,j))>=thres) then
                  write (fileout,'(i5,10f15.7)') i,(mat2(i,ind),ind=istart,iend)
                  exit
               endif
            enddo
         enddo
      else
      endif
   endif
end subroutine print_r2mat
