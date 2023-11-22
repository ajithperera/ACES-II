program get_acesmo
   implicit none
   integer fileinp,fileout
   character(len=80) keyword

   fileout = 6
   open (unit=fileout,file='alice.out',status='unknown')
   write (fileout,'(a,/)') '@ get_acesmo : Start program'

   call job_nwchem_mo(fileout)

   write (fileout,'(/,a)') '@ get_acesmo : End of program'
   close (fileout)
end program get_acesmo


subroutine job_nwchem_mo(fileout)
   use mod_movec
   use mod_info
   implicit none
   integer fileout,filemo
   integer nocc_sym(8),nvir_sym(8)
   integer ntot,nocc,nvir,nlin,ncart,nsph,ind1(1)
   integer natom,mxang
   integer uhf,ref,i,j 
   double precision, allocatable :: movec_ab(:,:,:)
   double precision, allocatable :: movec(:,:)
   double precision, allocatable :: trmat(:,:)
   character(len=3), allocatable :: nw_nuc(:)
   integer,          allocatable :: nw_ang(:,:)
   integer,          allocatable :: nw_ind(:,:)
   character(5) molab(8)
   integer nirrep,occupya(8),occupyb(8),numbasir(8),nlindep(8)
   logical yesno
   filemo=31

   !read ACES2 info (a2info.dat)
   call get_aces2_info(nirrep,occupya,occupyb,numbasir,nlindep,molab,fileout)

   !read info (nwchem.info)
   call info_read_mo(ntot,nocc,nvir,nlin,nirrep,nlindep,fileout,.true.)
   call info_symblk(nocc,nvir,nlin,nocc_sym,nvir_sym,nirrep,molab,fileout)
   call info_get_rotindex(nocc,nvir)

   !read 'nwchem.aolab'
   allocate(nw_nuc(ntot))
   allocate(nw_ind(ntot,3))
   allocate(nw_ang(ntot,3))
   call info_read_aolab(ntot,nw_nuc,nw_ind,nw_ang,fileout)

!  get NWChem MOs
   allocate(movec(ntot,ntot))
   allocate(movec_ab(ntot,ntot,2))
   call movec_read_mo(ntot,movec_ab,fileout,uhf)
   if (uhf==1) then
      call movec_write_mo_init(occupya,occupya,nirrep,filemo,fileout)
   else
      call movec_write_mo_init(occupya,occupyb,nirrep,filemo,fileout)
   endif

   do ref=1,uhf
      do i=1,ntot         
      do j=1,ntot         
         movec(i,j)=movec_ab(i,j,ref)
      enddo
      enddo

      !change MO ordering wrt symmetry block
      call movec_symblk(ntot,nocc,nvir,movec,rot_occ,rot_vir, &
           orbsym_occ,orbsym_vir,fileout)

      write(fileout,'(/,a,i1)') "* Print the orginal NWChem MO vector, spin=",ref 
      call print_mo(movec,ntot,ntot,ntot,nocc,nvir, &
           orbsym_occ,orbsym_vir,nw_nuc,aolab,5,fileout)

      write(fileout,'(/,a,i1)') "* Rotate MO vector, spin=",ref
      call job_rotate_ao(ntot,nocc,nvir,ref,movec,aolab, &
              orbsym_occ,orbsym_vir,nw_nuc,nw_ind,nw_ang,fileout)
      
      write(fileout,'(/,a,i1)') "* Write the new MO vector, spin=",ref
      call movec_write_mo_aces2(ntot,nocc,nvir,nocc_sym,nvir_sym, &
             movec,nirrep,numbasir,filemo,fileout)
   enddo 

   close(filemo)
   write(fileout,'(a)') '* New MO is written into NEWMOS.work'
   if (allocated(movec_ab)) deallocate(movec_ab)
   if (allocated(movec))    deallocate(movec)
end subroutine job_nwchem_mo



subroutine job_rotate_ao(ntot,nocc,nvir,ref,movec,aolab, &
              orbsym_occ,orbsym_vir,nw_nuc,nw_ind,nw_ang,fileout)
   use mod_movec
   implicit none
   integer ntot,nocc,nvir,ref,fileout
   double precision movec(ntot,ntot)
   double precision MatTr1(ntot,ntot)
   double precision MatTr2(ntot,ntot)
   double precision MatS(ntot,ntot)

   character(len=10) aolab(ntot)
   character(len=5) orbsym_occ(ntot),orbsym_vir(ntot)
   character(len=3) nw_nuc(ntot)
   integer nw_ind(ntot,3)
   integer nw_ang(ntot,3)

!  get_aces2_basisinfo
   integer aorot(ntot)
   double precision ao_norm(ntot)
   character(2) ao_nuc(ntot)
   character(5) ao_str(ntot)
   integer ao_ind(ntot,3),ao_ang(ntot,3),i
   logical loverlap_norm

!  get mo (from nwchem), overlap (from aces2) and calculate C'SC

!  ao_nuc/ao_str/ao_ind/ao_ang
   write(fileout,'(/,a)') 'Print read BASINFO'
   call get_aces2_basinfo(ntot,ao_nuc,ao_str,ao_ind,ao_ang, &
        'ACES2_BASINFO',fileout)

!  read S, AOSO, SOAO matrix from ACES2
   call get_aces2_matrix1(ntot,MatS,'a2info_overlap.dat',fileout)
   write(fileout,'(/,a)') '* Print the overlap matrix'
   call print_r2mat(MatS,ntot,ntot,5,fileout,1)
   write(fileout,'(/,a)') "* Overlap matrix (diag)"
   do i=1,ntot
      write(fileout,*) MatS(i,i)
   enddo 

   call get_aces2_matrix1(ntot,MatTr2,'a2info_mat_soao.dat',fileout)
   write(fileout,'(/,a)') '* Print the transformation matrix (SOAO)'
   call print_r2mat(MatTr2,ntot,ntot,5,fileout,2)

!  Rotato AO ordering (aorot(NWCHEM/ACES2) -> movec, aolab)
   call movec_match_aoind(ntot,ao_ind,ao_ang,ao_str,fileout,1)
   call movec_aorot_nw2ac(ntot,nw_ind,nw_ang,nw_nuc, &
        ao_ind,ao_ang,ao_nuc,aorot,fileout)

!  Rotate movec/aolab/nuc_name
   call movec_ao_rotate1(ntot,movec,aorot,aolab,nw_nuc,fileout,0)

!  Rotate AO ordering (SOAO)
!  call movec_ao_rotSOAO(ntot,movec,MatTr2,aolab,nw_nuc,fileout)

!  Normalize AO coefficient in (rotated) MO
   loverlap_norm = .false.
   if (loverlap_norm) then
      call movec_ao_norm(ntot,nw_ang,movec,MatS,aolab,aorot,ao_norm,fileout,1)
   else
      call movec_ao_norm(ntot,nw_ang,movec,MatS,aolab,aorot,ao_norm,fileout,2)
   endif

   write(fileout,'(/,a,i1)') "* Print AO rotated NWChem MO vector (step1), spin=",ref
   call print_mo(movec,ntot,ntot,ntot,nocc,nvir, &
        orbsym_occ,orbsym_vir,nw_nuc,aolab,5,fileout)

   call print_overlap(ntot,movec,movec,MatS,fileout,'CSC')
end subroutine job_rotate_ao 

