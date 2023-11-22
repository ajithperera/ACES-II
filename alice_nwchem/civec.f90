module mod_civec
   double precision, allocatable :: ntoval1(:)
   double precision, allocatable :: ntoval2(:)
   double precision, allocatable :: ntovec1(:,:)
   double precision, allocatable :: ntovec2(:,:)

   double precision, allocatable :: tmat1(:,:)
   double precision, allocatable :: tmat2(:,:)
contains



subroutine civec_tmat_rotate(nocc,nvir,rot_occ,rot_vir,fileout)
   use mod_print
   implicit none
   integer nocc,nvir,fileout
   integer rot_occ(nocc),rot_vir(nvir)
   integer i,a

   allocate(tmat2(nvir,nocc))
   do i = 1, nocc
      do a = 1, nvir
         tmat2(a,i)=tmat1(rot_vir(a),rot_occ(i))
      enddo
   enddo

   write(fileout,'(/,a)') "* Print Tai matrix (rotated w/sym.label)"
   call print_r2mat(tmat2,nvir,nocc,5,fileout,1) 
end subroutine civec_tmat_rotate


subroutine civec_read_tmat(ntot,nocc,nvir,nlin,orbeng1,orbeng2,fileout)
   use mod_print
   implicit none
   integer fileinp,fileout
   integer ntot,nocc,nvir,nlin,ind,i,a
   double precision orbeng1(nocc),orbeng2(nvir),thres(2)
   character(80) string

   thres(1) = -2.0d0
   thres(2) =  0.5d0

   fileinp = 11
   open(unit=fileinp,file='nwchem.civecs',status='old')
   do i = 1, 11
      read(fileinp,'(a)') string
   enddo

   allocate(tmat1(nvir,nocc))
   do i = 1, nocc
      do a = 1, nvir
         if (a <= nvir-nlin) then
            read(fileinp,*) tmat1(a,i)
         else
            tmat1(a,i) = 0.0d0
         endif
      enddo
   enddo

! subspace
   if (.true.) then
      write(fileout,'(a)') "civec_read_tmat: threshold"
      write(fileout,'(a,f15.7)') 'threshold,occ: ',thres(1)
      write(fileout,'(a,f15.7)') 'threshold,vir: ',thres(2)
      
      write(fileout,'(a)') "civec_read_tmat: orbital energy"
      do i = 1, nocc
         if (orbeng1(i).le.thres(1)) then
            write(fileout,'(i5,f15.7,a)') i,orbeng1(i),' *' 
         else
            write(fileout,'(i5,f15.7)')   i,orbeng1(i) 
         endif
      enddo
      do i = 1, nvir
         if (orbeng2(i).ge.thres(2)) then
            write(fileout,'(i5,f15.7,a)') i,orbeng2(i),' *' 
         else
            write(fileout,'(i5,f15.7)')   i,orbeng2(i) 
         endif
      enddo
      
      do i = 1, nocc
         do a = 1, nvir
            if (orbeng1(i).le.thres(1)) tmat1(a,i)=0.0d0 
            if (orbeng2(i).ge.thres(2)) tmat1(a,i)=0.0d0 
         enddo
      enddo
   endif
   close(fileinp)

   write(fileout,'(/,a)') "* Print Tai matrix (after reading from file)"
   call print_r2mat(tmat1,nvir,nocc,5,fileout,1) 
end subroutine civec_read_tmat


subroutine civec_get_transvec(nocc,nvir,nocc_sym,nvir_sym, &
           orbind2_occ,orbind2_vir,orbsym2_occ,orbsym2_vir,fileout)
   use mod_print
   implicit none
   integer nocc,nvir,nocc_sym(8),nvir_sym(8),fileout
   integer orbind2_occ(nocc),orbind2_vir(nvir)
   character(len=5) orbsym2_occ(nocc),orbsym2_vir(nvir)

   integer i,j,k,ind,ishift,filento
   logical lprint,yesno
   double precision Doo(nocc,nocc),Dvv(nvir,nvir),thres
   double precision, allocatable :: Mat1(:,:) 
   double precision, allocatable :: Mat2(:,:) 
   character(1) blank
   data blank /" "/

   allocate(ntoval1(nocc)) 
   allocate(ntoval2(nvir))
   allocate(ntovec1(nocc,nocc)) 
   allocate(ntovec2(nvir,nvir))

   call zero(ntoval1,nocc)
   call zero(ntoval2,nvir)
   call zero(ntovec1,nocc*nocc)
   call zero(ntovec2,nvir*nvir)

   lprint = .true.
   thres = 1.0D-7
   filento = 31

   inquire(file='NTO_EXCIT',exist=yesno)
   if(yesno) then
     open(unit=filento,file='NTO_EXCIT',status='old',form='formatted')
     close(unit=filento,status='delete')
   endif
   open(unit=filento,file='NTO_EXCIT',status='new',form='formatted')
  
   write(filento,'(a,x,i1)') "Symmetry of the excited State: ",0 
   write(filento,'(a)') blank

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
   ind = 0
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
            ntoval1(j+ishift) = Mat1(j,j)
            if (ntoval1(j+ishift)>=thres) then
               write(fileout,'(i5,f15.7)') j,ntoval1(j+ishift)
            endif
         enddo

         write(filento,'(a,x,i1)') "Symmetry of occ. orbitals:",i
         write(filento,'(a)') blank
         do j = 1, nocc_sym(i)
            ind = ind + 1
            write(filento,'(a,i7,x,f15.10)') "OCC",ind,Mat1(j,j)
         enddo
         write(filento,'(a)') blank

         do j = 1, nocc_sym(i)
            do k = 1, nocc_sym(i)
               ntovec1(j+ishift,k+ishift) = Mat2(j,k) 
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
   ind = 0
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
            ntoval2(j+ishift) = Mat1(j,j)
            if (ntoval2(j+ishift)>=thres) then
               write(fileout,'(i5,f15.7)') j,ntoval2(j+ishift)
            endif
         enddo

         write(filento,'(a,x,i1)') "Symmetry of vrt. orbitals:",i
         write(filento,'(a)') blank
         do j = 1, nvir_sym(i)
            ind = ind + 1
            write(filento,'(a,i7,x,f15.10)') "VRT",ind,Mat1(j,j)
         enddo
         write(filento,'(a)') blank


         do j = 1, nvir_sym(i)
            do k = 1, nvir_sym(i)
               ntovec2(j+ishift,k+ishift) = Mat2(j,k)
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
      call print_r2mat(ntoval1,nocc,1,5,fileout,2)
      write(fileout,'(/,a)') "Print Full NTO rotation vector (occupied)"
      call print_r2mat(ntovec1,nocc,nocc,5,fileout,2)
      write(fileout,'(/,a)') "Print Full NTO weight (virtual)"
      call print_r2mat(ntoval2,nvir,1,5,fileout,2)
      write(fileout,'(/,a)') "Print Full NTO rotation vector (virtual)"
      call print_r2mat(ntovec2,nvir,nvir,5,fileout,2)
   endif

   close(unit=filento,status='keep')
end subroutine civec_get_transvec


subroutine civec_mem_deallocation
   implicit none

   if (allocated(ntoval1)) deallocate(ntoval1)
   if (allocated(ntoval2)) deallocate(ntoval2)
   if (allocated(ntovec1)) deallocate(ntovec1)
   if (allocated(ntovec2)) deallocate(ntovec2)
   if (allocated(tmat1)) deallocate(tmat1)
   if (allocated(tmat2)) deallocate(tmat2)
end subroutine civec_mem_deallocation

end module mod_civec
