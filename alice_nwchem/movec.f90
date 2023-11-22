module mod_movec
   double precision, allocatable :: movec1(:,:)
   double precision, allocatable :: movec2(:,:)
!  double precision, allocatable :: nto(:,:)
contains

subroutine movec_read_mo(ntot,movec,fileout,uhf)
   use mod_print
   implicit none
   integer fileinp,fileout
   integer ntot,ntotao,ntotmo,nlin,nrest,nline1,i,j,k
   integer ref, uhf
   logical lprint
   character(80) string
   double precision  movec(ntot,ntot,2)
   double precision, allocatable :: mo_occup(:,:)
   double precision, allocatable :: mo_eng(:,:)

   fileinp = 21
   open(unit=fileinp,file='nwchem.movecs',status='old',position='rewind')
   do i = 1, 11
      read(fileinp,'(a80)') string
   enddo
   read(fileinp,*) uhf
   read(fileinp,*) ntotao !ao init.
   read(fileinp,*) ntotmo !ao without nlin
   nlin = ntotao-ntotmo
   nrest = mod(ntotao,3)
   nline1= (ntotao-nrest)/3
   

   allocate (mo_occup(ntotao,2))
   allocate (mo_eng(ntotao,2))

   do ref =1, uhf
   do i = 1, nline1
      read(fileinp,*) (mo_occup((i-1)*3+j,ref),j=1,3) 
   enddo
   if (nrest>=1) then
      read(fileinp,*) (mo_occup((nline1)*3+j,ref),j=1,nrest) 
   endif
   if (ref .eq.1) then
   write(fileout,'(/,a)') "* Print alpha mo occupancy"
   write(fileout,'(10(x,f5.2))') (mo_occup(i,ref),i=1,ntotao)
   else
   write(fileout,'(/,a)') "* Print beta  mo occupancy"
   write(fileout,'(10(x,f5.2))') (mo_occup(i,ref),i=1,ntotao)
   endif 

   do i = 1, nline1
      read(fileinp,*) (mo_eng((i-1)*3+j,ref),j=1,3)
   enddo
   if (nrest>=1) then
      read(fileinp,*) (mo_eng((nline1)*3+j,ref),j=1,nrest)
   endif

   do k = 1, ntotmo
      do i = 1, nline1
         read(fileinp,*) (movec((i-1)*3+j,k,ref),j=1,3)
      enddo
      if (nrest>=1) then
         read(fileinp,*) (movec((nline1)*3+j,k,ref),j=1,nrest)
      endif
   enddo
   if (nlin>=1) then
      do k = ntotmo+1,ntotao
         do j = 1, ntotao
            movec(j,k,ref) = 0.0d0
         enddo
      enddo
   endif
   enddo 

   do ref = 1, uhf
   if (ref.eq.1) then
   write(fileout,'(/,a)') "* movec_read_mo: Print alpha mo vectors"
   call print_r2mat(movec(1,1,ref),ntot,ntot,5,fileout,1)
   elseif (ref.eq.2) then
   write(fileout,'(/,a)') "* movec_read_mo: Print beta  mo vectors"
   call print_r2mat(movec(1,1,ref),ntot,ntot,5,fileout,1)
   endif
   enddo 

!  deallocate
   deallocate(mo_occup)
   deallocate(mo_eng)

   close(fileinp)
end subroutine movec_read_mo


subroutine movec_write_mo_aces2(ntot,nocc,nvir,nocc_sym,nvir_sym,mat2, &
              nirrep,occupya,numbasir,fileout)
   use mod_print
   implicit none
   integer ntot,nocc,nvir,nocc_sym(8),nvir_sym(8),fileout,fileread
   integer nirrep,occupya(8),numbasir(8)
   logical yesno
   integer irow,icol
   integer i,j,k,ind,rest,num2,iline,istart,iend
   double precision mat2(ntot,ntot)
   double precision,allocatable :: mat3(:,:)
   integer isym,dim0,ind1,ind2,ind3

   iline = 4
   fileread = 21
   inquire(file='NEWMOS.work',exist=yesno)
   if (yesno) then
      write (fileout,'(/,a)') ' NEWMOS.work already exists and will be deleted.'
      open(fileread,file='NEWMOS.work',status='old',access='sequential',form='formatted')
      close(fileread,status='delete')
   endif
   open(unit=fileread,file='NEWMOS.work',status='new',position='rewind')

   write(fileread,"(8(1x,i5))") (occupya(i), i=1, nirrep)
   write(fileread,"(8(1x,i5))") (occupya(i), i=1, nirrep)

   ind1 = 0
   ind2 = nocc
   ind3 = 0
   write(fileout,*) 'numbasir',numbasir
   do isym=1,nirrep
!     if (numbasir(isym)==0) cycle
      if (numbasir(isym)/=0) then
      dim0 = numbasir(isym)
      write(fileout,'(a,i5)')'numbasir=',dim0
      allocate(mat3(dim0,dim0))

!     make matrix for each symmetry      
      if (nocc_sym(isym)/=0) then
         write(fileout,'(a,i1,a,4(x,a,i5))') 'mo(occ,block=',isym,')', &
              'col',(ind1+1),' ~ ',(ind1+nocc_sym(isym)), &
              'row',(ind3+1),' ~ ',(ind3+dim0)
         do j = 1, nocc_sym(isym) 
            ind1 = ind1 + 1
            do i = 1, dim0 
               mat3(i,j) = mat2(i+ind3,ind1) 
            enddo
         enddo
      endif

      if (nvir_sym(isym)/=0) then
         write(fileout,'(a,i1,a,4(x,a,i5))') 'mo(vir,block=',isym,')', &
              'col',(ind2+1),' ~ ',(ind2+nvir_sym(isym)), &
              'row',(ind3+1),' ~ ',(ind3+dim0)
         do j = nocc_sym(isym)+1, nocc_sym(isym)+nvir_sym(isym) 
            ind2 = ind2 + 1
            do i = 1, dim0 
               mat3(i,j) = mat2(i+ind3,ind2) 
            enddo
         enddo
      endif

      ind3 = ind3+dim0

!     For temporary
      icol = nocc_sym(isym)+nvir_sym(isym) 
      irow = dim0

      write(fileout,'(a,2i5)') 'Start writing ',irow,icol
      
      ! icol = iline*num2 + rest
      rest = mod(icol,iline)
      num2 = (icol-rest)/iline
      do k = 1, num2
         istart = iline*(k-1) + 1
         iend = iline*(k-1) + iline
         do i = 1, irow
            write (fileread,'(4f20.10)') (mat3(i,ind),ind=istart,iend)
         end do
      end do
      
      if (rest>0) then
         istart = iline*num2+1
         iend = iline*num2+rest
         do i = 1, irow
            write (fileread,'(4f20.10)') (mat3(i,ind),ind=istart,iend)
         enddo
      endif

      deallocate(mat3)
      write(fileout,'(a,i1)') 'End symmetry,',isym
      endif
   enddo
   close(fileread)
   write(fileout,'(a)') 'NTO MO is written into NEWMOS.work'
end subroutine movec_write_mo_aces2


subroutine movec_symblk(ntot,nocc,nvir,movec1,rot_occ,rot_vir, &
           molab_occ,molab_vir,fileout)
   use mod_print
   implicit none
   integer ntot,nocc,nvir,fileout,i,j
   integer rot_occ(nocc),rot_vir(nvir)
   double precision movec1(ntot,ntot)
   double precision movec2(ntot,ntot)
   character(len=5) molab_occ(nocc)
   character(len=5) molab_vir(nvir)
   character(len=5) molab_tmp(ntot)

!  allocate(movec2(ntot,ntot))

!  occupied
   do i = 1, nocc
      do j = 1, ntot 
         movec2(j,i) = movec1(j,rot_occ(i))
      enddo
      molab_tmp(i) = molab_occ(rot_occ(i)) 
   enddo

!  virtual
   do i = 1, nvir
      do j = 1, ntot 
         movec2(j,nocc+i) = movec1(j,nocc+rot_vir(i))
      enddo
      molab_tmp(nocc+i) = molab_vir(rot_vir(i)) 
   enddo

!  copy back to movec1
   do i = 1, ntot
      do j = 1, ntot
         movec1(j,i) = movec2(j,i)
      enddo
      if (i>nocc) then
         molab_vir(i-nocc) = molab_tmp(i)
      else
         molab_occ(i) = molab_tmp(i)
      endif
   enddo
   
   write(fileout,'(/,a)') "*movec_symblk: Print mo vectors (symmetry blocked)"
   call print_r2mat(movec1,ntot,ntot,5,fileout,1)
end subroutine movec_symblk


subroutine movec_rotate_soao(ntot,Mat0,Mat1,Mat2,fileout)
   use mod_print
   implicit none
   integer ntot,fileout,i,j
   double precision Mat0(ntot,ntot),Mat1(ntot,ntot)
   double precision Mat2(ntot,ntot),Mat3(ntot,ntot)

!  call xgemm('T','N',ntot,ntot,ntot,1.0D0,Mat1,ntot,Mat0,ntot,0.0D0,Mat1,ntot)
!  call xgemm('N','N',ntot,ntot,ntot,1.0D0,Mat1,ntot,Mat1,ntot,0.0D0,Mat0,ntot)

!  call xgemm('N','T',ntot,ntot,ntot,1.0D0,Mat0,ntot,Mat2,ntot,0.0D0,Mat1,ntot)

!  call xgemm('T','N',ntot,ntot,ntot,1.0D0,Mat1,ntot,Mat2,ntot,0.0D0,Mat1,ntot)

!  call xgemm('N','N',ntot,ntot,ntot,1.0D0,Mat2,ntot,Mat1,ntot,0.0D0,Mat3,ntot)
!  write(fileout,'(/,a)') 'Print (SOAO)*(AOSO)'
!  call print_r2mat(Mat3,ntot,ntot,5,fileout,2)


   call xgemm('T','N',ntot,ntot,ntot,1.0D0,Mat2,ntot,Mat0,ntot,0.0D0,Mat3,ntot)
!  call xgemm('N','N',ntot,ntot,ntot,1.0D0,Mat3,ntot,Mat1,ntot,0.0D0,Mat0,ntot)

!  call xgemm('N','N',ntot,ntot,ntot,1.0D0,Mat2,ntot,Mat0,ntot,0.0D0,Mat1,ntot)
!  call xgemm('N','N',ntot,ntot,ntot,1.0D0,Mat1,ntot,Mat2,ntot,0.0D0,Mat0,ntot)

   if (.true.) then
      do i = 1, ntot
         do j = 1, ntot
            Mat0(i,j) = Mat3(i,j)
         enddo
      enddo
   endif

   write(fileout,'(/,a)') 'Print MO after transformation matrix (SOAO)'
   call print_r2mat(Mat0,ntot,ntot,5,fileout,2)
!  write(fileout,'(/,a)') 'Print Overlap after transformation matrix (SOAO)'
!  call print_r2mat(Mat1,ntot,ntot,5,fileout,1)
end subroutine movec_rotate_soao


subroutine movec_match_aoind(ntot,ao_ind,ao_ang,ao_str,fileout,bastype)
!  DATA (lable(I,1),I = 1,84) /'S000', 'P100', 'P010', 'P001',
!              'D200', 'D110', 'D101', 'D020', 'D011', 'D002',
!              'F300', 'F210', 'F201', 'F120', 'F111', 'F102',
!              'F030', 'F021', 'F012', 'F003',
!              'G400', 'G310', 'G301', 'G220', 'G211', 'G202',
!              'G130', 'G121', 'G112', 'G130', 'G040', 'G031',
!              'G022', 'G031', 'G004',
!              'H500', 'H410', 'H401', 'H320', 'H311', 'H302',
!              'H230', 'H221', 'H212', 'H203', 'H140', 'H131',
!              'H122', 'H113', 'H104', 'H050', 'H041', 'H032',
!              'H023', 'H014', 'H005',
!              'I600', 'I510', 'I501', 'I420', 'I411', 'I402',
!              'I330', 'I321', 'I312', 'I303', 'I240', 'I231',
!              'I222', 'I213', 'I204', 'I150', 'I141', 'I132',
!              'I123', 'I114', 'I105', 'I060', 'I051', 'I042',
!              'I033', 'I024', 'I015', 'I006'/
   use mod_fun
   implicit none
   integer ntot,fileout,bastype,i
   character(5) ao_str(ntot)
   integer ao_ind(ntot,3),ao_ang(ntot,3)
   do i = 1, ntot
      if (ao_ang(i,1)==1) then
         ao_ang(i,1) = 1
         call fun_ang_index1(ao_ang(i,2),1,ao_str(i))
!        call fun_ang_index2(ntot,i,ao_ind,ao_ang)
      elseif (ao_ang(i,1)==3) then
         ao_ang(i,1) = 2
         call fun_ang_index1(ao_ang(i,2),1,ao_str(i))
!        call fun_ang_index2(ntot,i,ao_ind,ao_ang)
      elseif (ao_ang(i,1)==6) then
         ao_ang(i,1) = 3
         call fun_ang_index1(ao_ang(i,2),2,ao_str(i))
!        call fun_ang_index2(ntot,i,ao_ind,ao_ang)
      elseif (ao_ang(i,1)==10) then
         ao_ang(i,1) = 4
         call fun_ang_index1(ao_ang(i,2),1,ao_str(i))
!        call fun_ang_index2(ntot,i,ao_ind,ao_ang)
      endif
   enddo
end subroutine movec_match_aoind


subroutine movec_get_aorot(ntot,nw_nuc,nw_ind,nw_ang,aorot, &
              ao_nuc,ao_str,ao_ind,ao_ang,fileout)
   implicit none
   integer ntot,fileout,i,j,itemp
   integer nw_ind(ntot,3),nw_ang(ntot,3)
   character(3) nw_nuc(ntot)
   character(2) ao_nuc(ntot)
   character(5) ao_str(ntot)
   integer ao_ind(ntot,3),ao_ang(ntot,3)
   integer aorot(ntot)

   write(fileout,'(/,a)') "* Print indices from NWCHEM and ACES2"
   do i = 1, ntot
!     write(fileout,'(a,4i4,a,a,5i3,x,a)') &
!          nw_nuc(i),nw_ind(i,1),(nw_ang(i,j),j=1,3), &
!          " : ",ao_nuc(i),(ao_ind(i,j),j=1,3),(ao_ang(i,j),j=1,2),ao_str(i)
      write(fileout,'(a,4i4,a,a,4i4,x,a)') &
           nw_nuc(i),nw_ind(i,1),(nw_ang(i,j),j=1,3), &
           " : ",ao_nuc(i),(ao_ind(i,j),j=1,1),(ao_ang(i,j),j=1,3),ao_str(i)
   enddo

   do i = 1, ntot
      aorot(i) = i
   enddo

   do i = 1, ntot-1
      do j = i+1, ntot
         if (nw_ind(i,1)==nw_ind(j,1)) then
            if (nw_ang(i,1)==nw_ang(j,1)) then
            if (nw_ang(i,2)< nw_ang(j,2)) then
               !switch i and j
               !write(fileout,'(a,2i5,2i10)')'switch ', &
               !     i,j,nw_ang(i,2),nw_ang(j,2)
               itemp = nw_ang(i,2)
               nw_ang(i,2) = nw_ang(j,2)
               nw_ang(j,2) = itemp
               itemp = aorot(i)
               aorot(i) = aorot(j)
               aorot(j) = itemp
               !write(fileout,'(a,2i5,2i10)')' ->    ', &
               !     aorot(i),aorot(j), &
               !     nw_ang(i,2),nw_ang(j,2)
            endif
            endif
         endif
      enddo
   enddo

   do i = 1, ntot-1
      do j = i+1, ntot
         if (nw_ind(i,1)==nw_ind(j,1)) then
            if (nw_ang(i,1)==nw_ang(j,1)) then
            if (nw_ang(i,2)==nw_ang(j,2)) then
               if (aorot(i)>aorot(j)) then
                  itemp = aorot(i)
                  aorot(i) = aorot(j)
                  aorot(j) = itemp
               endif
            endif
            endif
         endif
      enddo
   enddo

   write(fileout,'(a)') '* Print rotate vector between NWChem and ACES2'
   do i = 1, ntot
      write(fileout,'(4i5)') i,aorot(i), &
                             nw_ang(i,2),nw_ang(aorot(i),3)
   enddo
end subroutine movec_get_aorot


subroutine movec_rot_trmat(nsph,ncart,aoind2,trmat,fileout)
   use mod_fun
   use mod_print
   implicit none
   integer ncart,nsph,aoind2(ncart,6),fileout
   double precision trmat(nsph,ncart)
   double precision trmat2(nsph,ncart)
   integer i,j,i2

   do i = 1, ncart
      i2 = aoind2(i,6)
      do j = 1, nsph
         trmat2(j,i)=trmat(j,i2)
      enddo 
   enddo 

   do i = 1, ncart
      do j = 1, nsph
         trmat(j,i)=trmat2(j,i)
      enddo 
   enddo 

   write(fileout,'(/,a)') '* trmat (new)'
   call print_rmat(trmat,ncart*nsph,nsph,ncart,fileout,2)
end subroutine movec_rot_trmat


subroutine movec_aorot_nw2ac(ntot,nw_ind,nw_ang,nw_nuc, &
           ao_ind,ao_ang,ao_nuc,aorot,fileout)
   use mod_fun
   implicit none
   integer ntot,fileout,i,j,k
   integer nw_ind(ntot,3),nw_ang(ntot,3),itemp(3)
   integer tmp_ind(ntot,3),tmp_ang(ntot,3)
   character(3) nw_nuc(ntot),tmp_nuc(ntot)
   character(2) ao_nuc(ntot)
   integer ao_ind(ntot,3),ao_ang(ntot,3)
   integer aorot(ntot)

   write(fileout,'(/,a)') "* Print indices from NWCHEM and ACES2 (before)"
   do i = 1, ntot
!     write(fileout,'(a,4i4,a,a,5i3)') &
!          nw_nuc(i),nw_ind(i,1),(nw_ang(i,j),j=1,3), &
!          " : ",ao_nuc(i),(ao_ind(i,j),j=1,3),(ao_ang(i,j),j=1,2)
      write(fileout,'(a,4i4,a,i3,a,a,6i4)') &
           nw_nuc(i),nw_ind(i,1),(nw_ang(i,j),j=1,3), &
           " : ",i," : ",ao_nuc(i),(ao_ind(i,j),j=1,3),(ao_ang(i,j),j=1,3)
   enddo

   do i = 1, ntot
      aorot(i) = 0
      do j = 1, ntot
         if (nw_ind(j,1)==ao_ind(i,1)) then
         if (nw_ang(j,1)==ao_ang(i,1)) then
         if (nw_ang(j,2)==ao_ang(i,2)) then
         if (nw_ang(j,3)==ao_ang(i,3)) then
            nw_ind(j,2)=ao_ind(i,2)
            nw_ind(j,3)=ao_ind(i,3)

            aorot(i) = j
!           do k = 1, 3
!              tmp_ind(i,k) = nw_ind(j,k)
!              tmp_ang(i,k) = nw_ang(j,k)
!           enddo
!           tmp_nuc(i) = nw_nuc(j)
         endif
         endif
         endif
         endif
      enddo
   enddo

!  do i = 1, ntot
!     do k = 1, 3
!        nw_ind(i,k)=tmp_ind(i,k)
!        nw_ang(i,k)=tmp_ang(i,k)
!     enddo
!     nw_nuc(i)=tmp_nuc(i)
!  enddo

   write(fileout,'(/,a)') "* Print indices from NWCHEM and ACES2 (after)"
   do i = 1, ntot
      k = aorot(i)
      write(fileout,'(a,6i4,a,i3,a,a,6i4)') &
           nw_nuc(k),(nw_ind(k,j),j=1,3),(nw_ang(k,j),j=1,3), &
           " : ",k," : ",ao_nuc(i),(ao_ind(i,j),j=1,3),(ao_ang(i,j),j=1,3)
   enddo

end subroutine movec_aorot_nw2ac


subroutine movec_aorot_nw2ac2(ntot,nw_ind,nw_ang,nw_nuc, &
           ao_ind,ao_ang,ao_nuc,aorot,fileout)
   implicit none
   integer ntot,fileout,i,j,k
   integer nw_ind(ntot,3),nw_ang(ntot,3)
   character(3) nw_nuc(ntot)
   character(2) ao_nuc(ntot)
   integer ao_ind(ntot,3),ao_ang(ntot,3)
   integer aorot(ntot)


   do i = 1, ntot-1
      do j = i+1, ntot
         if (nw_ang(j,1)==ao_ang(i,1)) then
         if (nw_ang(j,2)==ao_ang(i,2)) then
         if (nw_ang(j,3)==ao_ang(i,3)) then
         if (nw_ind(j,2)==ao_ind(i,2)) then
!           aorot(i) = j
         endif
         endif
         endif
         endif
      enddo
   enddo

   write(fileout,'(/,a)') "* Print indices from NWCHEM and ACES2 (NW2AC2) "
   do i = 1, ntot
      k = aorot(i)
      write(fileout,'(a,6i4,a,i3,a,a,6i4)') &
           nw_nuc(i),(nw_ind(i,j),j=1,3),(nw_ang(i,j),j=1,3), &
           " : ",k," : ",ao_nuc(i),(ao_ind(i,j),j=1,3),(ao_ang(i,j),j=1,3)
   enddo
end subroutine movec_aorot_nw2ac2


subroutine movec_ao_rotate1(ntot,Mat1,aorot,aolab,nuc_name,fileout,itype)
! Rotate AO with aorot (matching NWChem orbital ordering to ACES2, for no-symmetry)
   implicit none
   integer ntot,aorot(ntot),fileout,itype,i,j
   double precision Mat1(ntot,ntot),Mat2(ntot,ntot)
   character(len=10) aolab(ntot),aolab2(ntot)
   character(len=3)  nuc_name(ntot),nuc_name2(ntot)

   do i = 1, ntot
      do j = 1, ntot
         Mat2(j,i) = Mat1(aorot(j),i)
      enddo
      aolab2(i) = aolab(aorot(i))
      nuc_name2(i) = nuc_name(aorot(i)) 
   enddo

   do i = 1, ntot
      do j = 1, ntot
         Mat1(j,i) = Mat2(j,i)
      enddo
      aolab(i) = aolab2(i)
      nuc_name(i) = nuc_name2(i)
   enddo
end subroutine movec_ao_rotate1


subroutine movec_ao_rotSOAO(ntot,movec,MatTr,aolab,nuc_name,fileout)
! Rotate AO with SOAO matrix (from ACES2, for symmetry)
   use mod_print
   implicit none
   integer ntot,aorot(ntot),fileout,i,j
   double precision movec(ntot,ntot),MatTr(ntot,ntot)
   double precision Mat1(ntot,ntot),Mat2(ntot,ntot)
   character(len=10) aolab(ntot),aolab2(ntot)
   character(len=3)  nuc_name(ntot),nuc_name2(ntot)
   double precision thres

   call xgemm('T','N',ntot,ntot,ntot,1.0D0,MatTr,ntot,MatTr,ntot,0.0D0,Mat2,ntot)
   write(fileout,'(/,a)') '* Print MatTr*MatTr'
   call print_r2mat(Mat2,ntot,ntot,5,fileout,2)
   

   thres = 0.1D-5
   write(fileout,'(/,a)') '* Print ao reordering (aolab)'
   do i = 1, ntot
      do j = 1, ntot
         if (abs(MatTr(j,i))>thres) then
            write(fileout,'(i5,a,i5,a10)') i,' -> ',j,aolab(j)
            aolab2(i) = aolab(j)
            nuc_name2(i) = nuc_name(j)
         endif
      enddo
   enddo

   do i = 1, ntot
      aolab(i) = aolab2(i)
      nuc_name(i) = nuc_name2(i)
   enddo
   
   call xgemm('T','N',ntot,ntot,ntot,1.0D0,MatTr,ntot,movec,ntot,0.0D0,Mat1,ntot)

!  call dcopy(ntot*ntot,Mat1,1,movec,1)
   do i = 1, ntot
      do j = 1, ntot
!        movec(i,j) = Mat1(i,j)
         movec(i,j) = Mat1(i,j)/Mat2(i,i)
      enddo
   enddo
end subroutine movec_ao_rotSOAO


subroutine movec_ao_norm1(ntot,nw_ang,ao_norm,fileout)
   implicit none
   integer ntot,fileout
   integer nw_ang(ntot,3)
   double precision ao_norm(ntot)
   integer i,itemp

   do i = 1, ntot
      ao_norm(i) = 1.0D0
      itemp = nw_ang(i,2)
      if (nw_ang(i,1)==3) then
         if (itemp==110 .or. itemp==101 .or. itemp==11) then
            ao_norm(i) = dsqrt(1.0D0/3.0D0)
!ycp-b
         elseif (itemp==200 .or. itemp==20 .or. itemp==2) then
            ao_norm(i) = dsqrt(1.0D0/3.0D0)
!           ao_norm(i) = 1.0D0/3.0D0
!ycp-e
         endif
      elseif (nw_ang(i,1)==4) then
         if (itemp==210 .or. itemp==201 .or. itemp==120 .or. &
             itemp==102 .or. itemp==21  .or. itemp==12) then
!           ao_norm(i) = dsqrt(1.0D0/5.0D0)/3.0D0
            ao_norm(i) = dsqrt(1.0D0/15.0D0)
         elseif (itemp==111) then
            ao_norm(i) = dsqrt(1.0D0/15.0D0)
         elseif (itemp==300 .or. itemp==30 .or. itemp==3) then 
!           ao_norm(i) = 1.0D0/15.0D0
            ao_norm(i) = dsqrt(1.0D0/15.0D0)
         endif
      endif
   enddo
end subroutine movec_ao_norm1


subroutine movec_ao_norm(ntot,Mat1,MatS,aolab,aorot,ao_norm, &
                           fileout,itype) 
   implicit none
   integer ntot,aorot(ntot),fileout,itype,i,j
   character(len=10) aolab(ntot)
   double precision Mat1(ntot,ntot),Mat2(ntot,ntot),MatS(ntot,ntot)
   double precision ao_norm(ntot)

   if (itype==1) then
      !correct with overlap matrix
      do i = 1, ntot
         ao_norm(aorot(i)) = ao_norm(aorot(i))/dsqrt(MatS(i,i))
!        ao_norm(aorot(i)) = ao_norm(aorot(i))
         write(fileout,*) 'MatS',i,MatS(i,i)
      enddo
   endif

   write(fileout,'(/,a)') '* Print normalization coefficients'
   do i = 1, ntot
      write(fileout,'(i5,x,a5,f15.7)') i,aolab(i),ao_norm(aorot(i))
   enddo

   do i = 1, ntot
      do j = 1, ntot
         Mat1(j,i) = ao_norm(aorot(j)) * Mat1(j,i)
      enddo
   enddo
end subroutine movec_ao_norm


subroutine movec_get_nto(ntot,nocc,nvir,ntoval1,ntoval2,ntovec1,ntovec2, &
                         movec,nto,fileout)
   use mod_print
   implicit none
   integer ntot,nocc,nvir,fileout
   double precision ntoval1(nocc),ntovec1(nocc,nocc)
   double precision ntoval2(nvir),ntovec2(nvir,nvir)
   double precision movec(ntot,ntot),nto(ntot,ntot)
   integer i,j,k
   double precision temp

   call zero(nto,ntot*ntot)

   !occupied, nto = movec * ntovec1 
   do i = 1, ntot   
      do j = 1, nocc
         temp = 0.0d0
         do k = 1, nocc
            temp = temp + movec(i,k) * ntovec1(k,j) 
         enddo
         nto(i,j) = temp
      enddo
   enddo

   !virtual, nto = movec * ntovec1
   do i = 1, ntot
      do j = 1, nvir
         temp = 0.0d0
         do k = 1, nvir
            temp = temp + movec(i,nocc+k) * ntovec2(k,j)
         enddo
         nto(i,nocc+j) = temp
      enddo
   enddo

   write(fileout,'(/,a)') "* Print NTO vectors (symmetry blocked)"
   call print_r2mat(nto,ntot,ntot,5,fileout,1)
end subroutine movec_get_nto


subroutine movec_mem_deallocation
   implicit none
   if (allocated(movec1)) deallocate(movec1)
   if (allocated(movec2)) deallocate(movec2)
!  if (allocated(nto)) deallocate(nto)
end subroutine movec_mem_deallocation

end module mod_movec
