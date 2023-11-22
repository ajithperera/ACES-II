subroutine civec_get_transvec2(tmat2,nocc,nvir,nocc_sym,nvir_sym, &
           orbind2_occ,orbind2_vir,orbsym2_occ,orbsym2_vir,fileout)
   use mod_print
   implicit none
   integer nocc,nvir,nocc_sym(8),nvir_sym(8),fileout
   integer orbind2_occ(nocc),orbind2_vir(nvir)
   character(len=5) orbsym2_occ(nocc),orbsym2_vir(nvir)

   integer i,j,k,ishift
   logical lprint
   double precision Doo(nocc,nocc),Dvv(nvir,nvir),tmat2(nvir,nocc),thres
   double precision tmat3(nvir,nocc)
   double precision, allocatable :: Mat1(:,:)
   double precision, allocatable :: Mat2(:,:)
   double precision, allocatable :: eigval1(:)
   double precision, allocatable :: eigval2(:)
   double precision, allocatable :: eigvec1(:,:)
   double precision, allocatable :: eigvec2(:,:)

   write(fileout,'(/,a)') "* Print number of each symmetry (occupied)"
   write(fileout,'(8i5)') (nocc_sym(i),i=1,8)
   write(fileout,'(/,a)') "* Print symmetry lable (occupied)"
   write(fileout,'(5(i5,a,a5))') (orbind2_occ(i),":",orbsym2_occ(i),i=1,nocc)
   write(fileout,'(/,a)') "* Print number of each symmetry (virtual)"
   write(fileout,'(8i5)') (nvir_sym(i),i=1,8)
   write(fileout,'(/,a)') "* Print symmetry lable (virtual)"
   write(fileout,'(5(i5,a,a5))') (orbind2_vir(i),":",orbsym2_vir(i),i=1,nvir)

   write(fileout,'(/,a)') 'Print Tai (before rotate)'
   call print_r2mat(tmat2,nvir,nocc,5,fileout,1)

   do i = 1, nvir
      do j = 1, nocc
         tmat3(i,j)=tmat2(orbind2_vir(i),orbind2_occ(j))
      enddo
   enddo

   write(fileout,'(/,a)') 'Print Tai (after rotate)'
   call print_r2mat(tmat3,nvir,nocc,5,fileout,1)

!  For Doo
   call dgemm('T','N',nocc,nocc,nvir,1.0D0,tmat2,nvir,tmat2,nvir, &
              0.0D0,Doo,nocc)

   allocate(eigval1(nocc))
   allocate(eigval2(nvir))
   allocate(eigvec1(nocc,nocc))
   allocate(eigvec2(nvir,nvir))

   call zero(eigval1,nocc)
   call zero(eigval2,nvir)
   call zero(eigvec1,nocc*nocc)
   call zero(eigvec2,nvir*nvir)

   lprint = .true.
   thres = 1.0D-7

   ishift=0
   do i = 1, 8
      if (nocc_sym(i)>=1) then
         allocate(Mat1(nocc_sym(i),nocc_sym(i)))
         allocate(Mat2(nocc_sym(i),nocc_sym(i)))
         do j = 1, nocc_sym(i)
            do k = 1, nocc_sym(i)
               Mat1(j,k) = Doo(j+ishift,k+ishift)
            enddo
         enddo
         write(fileout,'(/,a,i)') "Print Mat1 (before diag.), sym=",i
         call print_r2mat(Mat1,nocc_sym(i),nocc_sym(i),5,fileout,1)

         call eig(Mat1,Mat2,1,nocc_sym(i),-1)
         write(fileout,'(/,a,i)') "Print Mat1: eigenvalue, sym=",i
         !call print_r2mat(Mat1,nocc_sym(i),nocc_sym(i),5,fileout,1)
         do j = 1, nocc_sym(i)
            eigval1(j+ishift) = Mat1(j,j)
            if (eigval1(j+ishift)>=thres) then
               write(fileout,'(i5,f15.7)') j,eigval1(j+ishift)
            endif
         enddo

         do j = 1, nocc_sym(i)
            do k = 1, nocc_sym(i)
               eigvec1(j+ishift,k+ishift) = Mat2(j,k)
            enddo
         enddo
         write(fileout,'(/,a,i)') "Print Mat2: eigenvector, sym=",i
         call print_r2mat(Mat2,nocc_sym(i),nocc_sym(i),5,fileout,1)

         deallocate(Mat1)
         deallocate(Mat2)
      endif
      ishift = ishift + nocc_sym(i)
   enddo

   deallocate(eigval1)
   deallocate(eigval2)
   deallocate(eigvec1)
   deallocate(eigvec2)

   return

! ---------------
   allocate(eigval1(nocc))
   allocate(eigval2(nvir))
   allocate(eigvec1(nocc,nocc))
   allocate(eigvec2(nvir,nvir))

   call zero(eigval1,nocc)
   call zero(eigval2,nvir)
   call zero(eigvec1,nocc*nocc)
   call zero(eigvec2,nvir*nvir)

   lprint = .true.
   thres = 1.0D-7

!  For Doo
   call dgemm('T','N',nocc,nocc,nvir,1.0D0,tmat2,nvir,tmat2,nvir, &
              0.0D0,Doo,nocc)

   if (lprint) then
      write(fileout,'(/,a)') "* Print number of each symmetry (occupied)"
      write(fileout,'(8i5)') (nocc_sym(i),i=1,8)
      write(fileout,'(/,a)') "* Print symmetry lable (occupied)"
      write(fileout,'(5(i5,a,a5))') (orbind2_occ(i),":",orbsym2_occ(i),i=1,nocc)
      write(fileout,'(/,a)') 'Print Doo'
      call print_r2mat(Doo,nocc,nocc,5,fileout,1)
   endif

   ishift=0
   do i = 1, 8
      if (nocc_sym(i)>=1) then
         allocate(Mat1(nocc_sym(i),nocc_sym(i)))
         allocate(Mat2(nocc_sym(i),nocc_sym(i)))
         do j = 1, nocc_sym(i)
            do k = 1, nocc_sym(i)
               Mat1(j,k) = Doo(j+ishift,k+ishift)
            enddo
         enddo
         call eig(Mat1,Mat2,1,nocc_sym(i),-1)
         write(fileout,'(/,a,i)') "Print Mat1: eigenvalue, sym=",i
         !call print_r2mat(Mat1,nocc_sym(i),nocc_sym(i),5,fileout,1)
         do j = 1, nocc_sym(i)
            eigval1(j+ishift) = Mat1(j,j)
            if (eigval1(j+ishift)>=thres) then
               write(fileout,'(i5,f15.7)') j,eigval1(j+ishift)
            endif
         enddo

         do j = 1, nocc_sym(i)
            do k = 1, nocc_sym(i)
               eigvec1(j+ishift,k+ishift) = Mat2(j,k)
            enddo
         enddo
         write(fileout,'(/,a,i)') "Print Mat2: eigenvector, sym=",i
         call print_r2mat(Mat2,nocc_sym(i),nocc_sym(i),5,fileout,1)

         deallocate(Mat1)
         deallocate(Mat2)
      endif
      ishift = ishift + nocc_sym(i)
   enddo

!  For Dvv
   call dgemm('N','T',nvir,nvir,nocc,1.0D0,tmat2,nvir,tmat2,nvir, &
              0.0D0,Dvv,nvir)

   if (lprint) then
      write(fileout,'(/,a)') "* Print number of each symmetry (virtual)"
      write(fileout,'(8i5)') (nvir_sym(i),i=1,8)
      write(fileout,'(/,a)') "* Print symmetry lable (virtual)"
      write(fileout,'(5(i5,a,a5))') (orbind2_vir(i),":",orbsym2_vir(i),i=1,nvir)
      write(fileout,'(/,a)') 'Print Dvv'
      call print_r2mat(Dvv,nvir,nvir,5,fileout,1)
   endif

   ishift=0
   do i = 1, 8
      if (nvir_sym(i)>=1) then
         allocate(Mat1(nvir_sym(i),nvir_sym(i)))
         allocate(Mat2(nvir_sym(i),nvir_sym(i)))
         do j = 1, nvir_sym(i)
            do k = 1, nvir_sym(i)
               Mat1(j,k) = Dvv(j+ishift,k+ishift)
            enddo
         enddo
         call eig(Mat1,Mat2,1,nvir_sym(i),-1)
         write(fileout,'(/,a,i)') "Print Mat1: eigenvalue, sym=",i
         !call print_r2mat(Mat1,nvir_sym(i),nvir_sym(i),5,fileout,1)
         do j = 1, nvir_sym(i)
            eigval2(j+ishift) = Mat1(j,j)
            if (eigval2(j+ishift)>=thres) then
               write(fileout,'(i5,f15.7)') j,eigval2(j+ishift)
            endif
         enddo

         do j = 1, nvir_sym(i)
            do k = 1, nvir_sym(i)
               eigvec2(j+ishift,k+ishift) = Mat2(j,k)
            enddo
         enddo

         write(fileout,'(/,a,i)') "Print Mat2: eigenvector, sym=",i
         call print_r2mat(Mat2,nvir_sym(i),nvir_sym(i),5,fileout,1)
         deallocate(Mat1)
         deallocate(Mat2)
      endif
      ishift = ishift + nvir_sym(i)
   enddo

   if (lprint) then
      write(fileout,'(/,a)') "Print Full NTO weight (occupied)"
      call print_r2mat(eigval1,nocc,1,5,fileout,2)
      write(fileout,'(/,a)') "Print Full NTO rotation vector (occupied)"
      call print_r2mat(eigvec1,nocc,nocc,5,fileout,2)
      write(fileout,'(/,a)') "Print Full NTO weight (virtual)"
      call print_r2mat(eigval2,nvir,1,5,fileout,2)
      write(fileout,'(/,a)') "Print Full NTO rotation vector (virtual)"
      call print_r2mat(eigvec2,nvir,nvir,5,fileout,2)
   endif
   
end subroutine civec_get_transvec2
