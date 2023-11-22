module mod_print
contains

subroutine print_rmat(mat,dim0,irow,icol,fileout,itype)
   implicit none
   integer dim0,irow,icol,fileout,itype
   integer i,j,k,ind,rest,num2,iline,istart,iend
   double precision mat(dim0)
   double precision, allocatable :: mat2(:,:)

   iline = 5

   if (itype==1) then
      do i = 1, irow
         write (fileout,'(i5,f15.7)') i, mat(i)
      end do
   else if (itype==2) then
      allocate (mat2(irow,icol))
      ind = 0
      do j = 1, icol
         do i = 1, irow
            ind = ind + 1
            mat2(i,j) = mat(ind)
         end do
      end do

      call print_r2mat(mat,irow,icol,iline,fileout,2)
      deallocate (mat2)
   end if
end subroutine print_rmat

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
               if (mat2(i,j)>=thres) then
                  write (fileout,'(i5,10f15.7)') i,(mat2(i,ind),ind=istart,iend)
                  exit
               endif
            enddo
         enddo
      else
      endif
   endif
end subroutine print_r2mat



subroutine print_rmat_symlab(mat,dim0,irow,icol,nocc,nvir, &
           symlab_occ,symlab_vir,fileout,itype)
   implicit none
   integer dim0,irow,icol,nocc,nvir,fileout,itype
   integer i,j,k,ind,rest,num2,iline,istart,iend
   character(len=5) symlab_occ(nocc),symlab_vir(nvir)
   double precision mat(dim0)
   double precision, allocatable :: mat2(:,:)

   iline = 5

   if (itype==1) then
      do i = 1, irow
         write (fileout,'(i5,f15.7)') i, mat(i)
      end do
   else if (itype==2) then
      allocate (mat2(irow,icol))
      ind = 0
      do i = 1, irow
         do j = 1, icol
            ind = ind + 1
            mat2(i,j) = mat(ind)
         end do
      end do
      call print_r2mat_symlab(mat,irow,icol,nocc,nvir, &
           symlab_occ,symlab_vir,iline,fileout)
      deallocate (mat2)
   endif
end subroutine print_rmat_symlab


subroutine print_r2mat_symlab(mat2,irow,icol,nocc,nvir, &
           symlab_occ,symlab_vir,iline,fileout)
   implicit none
   integer irow,icol,nocc,nvir,iline,fileout
   integer i,j,k,ind,rest,num2,istart,iend
   character(len=5) symlab_occ(nocc),symlab_vir(nvir)
   double precision mat2(irow,icol)

   ! icol = iline*num2 + rest
   rest = mod(icol,iline)
   num2 = (icol-rest)/iline
   do k = 1, num2
      istart = iline*(k-1) + 1
      iend = iline*(k-1) + iline
!     write (fileout,'(/,10x,10(4x,i5,a,5a))') & 
!           (i,':',symlab_occ(i),i=istart,iend)
      write (fileout,'(/,10x,10(5x,i4,a,a5))') & 
            (i,':',symlab_occ(i),i=istart,iend)
      do i = 1, irow
         write (fileout,'(i4,a,a5,10f15.7)') i,':',symlab_vir(i), &
               (mat2(i,ind),ind=istart,iend)
      end do
   end do

   if (rest>0) then
      istart = iline*num2+1
      iend = iline*num2+rest
      write (fileout,'(/,10x,10(4x,i4,a,a5))') & 
            (ind,':',symlab_occ(ind),ind=istart,iend)
      do i = 1, irow
         write (fileout,'(i4,a,a5,10f15.7)') i,':',symlab_vir(i), &
               (mat2(i,ind),ind=istart,iend)
      end do
   end if
end subroutine print_r2mat_symlab


subroutine print_mo(mat2,irow,icol,ntot,nocc,nvir, &
           symlab_occ,symlab_vir,nuc_name,aolab,iline,fileout)
   implicit none
   integer irow,icol,ntot,nocc,nvir,iline,fileout
   integer i,j,k,ind,rest,num2,istart,iend
   character(len=*) symlab_occ(nocc)
   character(len=*) symlab_vir(nvir)
   character(len=5) symlab(ntot)
   character(len=3) nuc_name(ntot)
   character(len=*) aolab(ntot)
   double precision mat2(irow,icol)

   ind = 0
   do i = 1, nocc 
      ind = ind + 1
      symlab(ind) = symlab_occ(i) 
   enddo
   do i = 1, nvir 
      ind = ind + 1
      symlab(ind) = symlab_vir(i) 
   enddo

   ! icol = iline*num2 + rest
   rest = mod(icol,iline)
   num2 = (icol-rest)/iline
   do k = 1, num2
      istart = iline*(k-1) + 1
      iend = iline*(k-1) + iline
      write (fileout,'(/,16x,10(2x,i4,a,a5))') &
            (ind,':',symlab(ind),ind=istart,iend)
      do i = 1, irow
         write (fileout,'(i4,x,a3,a7,10f12.7)') i,nuc_name(i),aolab(i), &
               (mat2(i,ind),ind=istart,iend)
      end do
   end do

   if (rest>0) then
      istart = iline*num2+1
      iend = iline*num2+rest
      write (fileout,'(/,16x,10(2x,i4,a,a5))') &
            (ind,':',symlab(ind),ind=istart,iend)
      do i = 1, irow
         write (fileout,'(i4,x,a3,a7,10f12.7)') i,nuc_name(i),aolab(i), &
               (mat2(i,ind),ind=istart,iend)
      end do
   end if
end subroutine print_mo

subroutine print_overlap(NAO,CMO1,CMO2,MatS,fileout,strovl)
   implicit none
   integer NAO,fileout
   character (len=*) strovl
   double precision :: CMO1(NAO,NAO)
   double precision :: CMO2(NAO,NAO)
   double precision :: MatS(NAO,NAO)

   double precision :: Mat1(NAO,NAO)
   double precision :: Mat2(NAO,NAO)
   integer i

   if (strovl=='C') then
      write (fileout,'(/,a)') "* Print |C1>"
      call print_rmat(CMO1,NAO*NAO,NAO,NAO,fileout,2)
      write (fileout,'(/,a)') "* Print |C2>"
      call print_rmat(CMO2,NAO*NAO,NAO,NAO,fileout,2)

   elseif (strovl=='CC') then
      write (fileout,'(/,a)') "* Print <C1|C2>"
      call zero(Mat1,NAO*NAO)
      call xgemm('t','n',NAO,NAO,NAO,1.0D0,CMO1,NAO,CMO2,NAO,0.0D0,Mat1,NAO)
      call print_rmat(Mat1,NAO*NAO,NAO,NAO,fileout,2)

   elseif (strovl=='CSC') then
      write (fileout,'(/,a)') "* Print <C1|S|C2>"
      call zero(Mat1,NAO*NAO)
      call zero(Mat2,NAO*NAO)
      call xgemm('T','N',NAO,NAO,NAO,1.0D0,CMO1,NAO,MatS,NAO,0.0D0,Mat1,NAO)
      call xgemm('N','N',NAO,NAO,NAO,1.0D0,Mat1,NAO,CMO2,NAO,0.0D0,Mat2,NAO)
      call print_rmat(Mat2,NAO*NAO,NAO,NAO,fileout,2)
      write(fileout,'(/,a)') "* Print diag."
      do i = 1, nao
         write(fileout,'(i5,f15.7,f20.15)') i,Mat2(i,i),Mat2(i,i)
      enddo

   endif
end subroutine print_overlap


end module mod_print
