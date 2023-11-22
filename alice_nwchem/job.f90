module mod_job
   double precision, allocatable :: MatS(:,:)
contains

subroutine job_control(keyword,fileout)
   implicit none
   integer fileout
   character(len=*) keyword

   if (keyword.eq."NWCHEM.COMPARE.ACES2.AO" .or. &
       keyword.eq."NWCHEM.NTO.CARTESIAN") then
      write(fileout,'(a,a)') "* Keyword = ",keyword
      call job_nwchem_nto(fileout,1) 
   elseif (keyword.eq."NWCHEM.NTO.SPHERICAL") then
      write(fileout,'(a,a)') "* Keyword = ",keyword
      call job_nwchem_nto(fileout,2) 
!  elseif (keyword.eq."NWCHEM.NTO") then
!     write(fileout,'(a,a)') "* Keyword = ",keyword
!     call nwchem_civec(fileout)
   elseif (keyword.eq."NWCHEM.MO.CARTESIAN") then
      write(fileout,'(a,a)') "* Keyword = ",keyword
      call job_nwchem_mo(fileout,1) 
   endif
end subroutine job_control


subroutine job_nwchem_mo_transform(natom,mxang,ncart,nsph,trmat,fileout)
   use mod_print
   use mod_aovec
   use mod_movec
   use mod_info
   use mod_get_aces2
   implicit none
   integer ncart,nsph,fileout
   integer natom,mxang
   double precision :: trmat(nsph,nsph)
   double precision :: trmat1(nsph,nsph)
   double precision :: trmat2(nsph,nsph)
   character(4), allocatable :: aostr(:,:)
   integer, allocatable :: aoind(:,:)
   integer, allocatable :: aoind2(:,:)
   integer, allocatable :: geoind(:,:)

   allocate(aostr(nsph,2))
   allocate(aoind(nsph,5))
   allocate(aoind2(ncart,6))
   allocate(geoind(natom,3))

   call info_read_acesgeo(natom,geoind,fileout)
   call info_read_aolab3(nsph,natom,mxang,aoind,aostr,fileout)

   call aovec_aosph_nw2ac(nsph,natom,aoind,aostr,geoind,trmat1,fileout)
   call aovec_trans_sym(nsph,trmat2,fileout)
   write(fileout,'(a)') 'test'
!  call aovec_trans_total(nsph,nsph,nsph,trmat1,trmat2,trmat)

!  call aovec_sph2cart_addfn(ncart,nsph,aoind,aoind2,fileout)
!  call aovec_make_trmat(nsph,ncart,aoind2,trmat,fileout)
!  call movec_rot_trmat(nsph,ncart,aoind2,trmat,fileout)

   deallocate(geoind)
   deallocate(aostr)
   deallocate(aoind)
   deallocate(aoind2)
end subroutine job_nwchem_mo_transform


subroutine job_nwchem_mo(fileout,bastype)
   use mod_print
   use mod_info
   use mod_aovec
   use mod_movec
   use mod_civec
   use mod_get_aces2
   implicit none
   integer fileout,bastype
   integer nocc_sym(8),nvir_sym(8)
   integer ntot,nocc,nvir,nlin,ncart,nsph,ind1(1)
   integer natom,mxang
   integer uhf,ref 
   double precision, allocatable :: movec(:,:,:)
   double precision, allocatable :: trmat(:,:)
   character(len=3), allocatable :: nw_nuc(:)
   integer,          allocatable :: nw_ang(:,:)
   integer,          allocatable :: nw_ind(:,:)
   character(5) molab(8)
   integer nirrep,occupya(8),numbasir(8),nlindep(8)
   logical yesno

   if (bastype==1) then
      !read ACES2 info (alice_aces2_info.txt)
      call get_aces2_info(nirrep,occupya,numbasir,nlindep,molab,fileout)

      !read info (nwchem.info)
      call info_read_mo(ntot,nocc,nvir,nlin,nirrep,nlindep,fileout,.true.)
      call info_symblk(nocc,nvir,nlin,nocc_sym,nvir_sym,nirrep,molab,fileout)
      call info_get_rotindex(nocc,nvir)

      !get mo
      allocate(movec(ntot,ntot,2))
      call movec_read_mo(ntot,movec,fileout,uhf)

      !change MO ordering wrt symmetry block
      do ref=1,uhf
      call movec_symblk(ntot,nocc,nvir,movec(1,1,ref),rot_occ,rot_vir, &
           orbsym_occ,orbsym_vir,fileout)
      enddo

      !read 'nwchem.aolab'
      allocate(nw_nuc(ntot))
      allocate(nw_ind(ntot,3))
      allocate(nw_ang(ntot,3))
      call info_read_aolab(ntot,nw_nuc,nw_ind,nw_ang,fileout)
      if (.true.) then
         ! MO (NWChem)
         do ref =1,uhf
         if (ref .eq. 1) then
         write(fileout,'(/,a,x,a)') "* Print alpha MO vector", &
              "(Original NWChem ordering w/AO lable)"
         call print_mo(movec(1,1,ref),ntot,ntot,ntot,nocc,nvir, &
              orbsym_occ,orbsym_vir,nw_nuc,aolab,5,fileout)
         elseif (ref .eq. 2) then 
         write(fileout,'(/,a,x,a)') "* Print beta MO vector", &
              "(Original NWChem ordering w/AO lable)"
         call print_mo(movec(1,1,ref),ntot,ntot,ntot,nocc,nvir, &
              orbsym_occ,orbsym_vir,nw_nuc,aolab,5,fileout)
         endif 
         enddo 
      endif

      !rotate AO
      do ref=1,uhf
      call job_rotate_ao(ntot,nocc,nvir,movec(1,1,ref),aolab, &
              orbsym_occ,orbsym_vir,nw_nuc,nw_ind,nw_ang,fileout)

      !nocc_sym(8), nvir_sym(8)
      call movec_write_mo_aces2(ntot,nocc,nvir,nocc_sym,nvir_sym, &
             movec(1,1,ref),nirrep,occupya,numbasir,fileout)
      enddo 
   endif
end subroutine job_nwchem_mo



subroutine job_nwchem_nto(fileout,bastype)
   use mod_print
   use mod_info
   use mod_aovec
   use mod_movec
   use mod_civec
   use mod_get_aces2
   implicit none
   integer fileout,bastype
   integer nocc_sym(8),nvir_sym(8)
   integer ntot,nocc,nvir,nlin,ncart,nsph,ind1(1)
   integer natom,mxang,uhf
   double precision, allocatable :: movec(:,:)
   double precision, allocatable :: trmat(:,:)
   character(len=3), allocatable :: nw_nuc(:)
   integer,          allocatable :: nw_ang(:,:)
   integer,          allocatable :: nw_ind(:,:)
   character(5) molab(8)
   integer nirrep,occupya(8),numbasir(8),nlindep(8)
   logical yesno

   if (bastype==1) then
      !read ACES2 info (alice_aces2_info.txt)
      call get_aces2_info(nirrep,occupya,numbasir,nlindep,molab,fileout)
      
      !read info (nwchem.info)
      call info_read_mo(ntot,nocc,nvir,nlin,nirrep,nlindep,fileout,.true.)
      call info_symblk(nocc,nvir,nlin,nocc_sym,nvir_sym,nirrep,molab,fileout)
      call info_get_rotindex(nocc,nvir)
      
      !get mo
      allocate(movec(ntot,ntot))
      call movec_read_mo(ntot,movec,fileout,uhf)
      
      !change MO ordering wrt symmetry block
      call movec_symblk(ntot,nocc,nvir,movec,rot_occ,rot_vir, &
           orbsym_occ,orbsym_vir,fileout)
      
      !read 'nwchem.aolab'
      allocate(nw_nuc(ntot))
      allocate(nw_ind(ntot,3))
      allocate(nw_ang(ntot,3))
      call info_read_aolab(ntot,nw_nuc,nw_ind,nw_ang,fileout)
      if (.true.) then
         ! MO (NWChem) 
         write(fileout,'(/,a,x,a)') "* Print MO vector", &
              "(Original NWChem ordering w/AO lable)"
         call print_mo(movec,ntot,ntot,ntot,nocc,nvir, &
              orbsym_occ,orbsym_vir,nw_nuc,aolab,5,fileout)
      endif
      
      !rotate AO
      call job_rotate_ao(ntot,nocc,nvir,movec,aolab, &
              orbsym_occ,orbsym_vir,nw_nuc,nw_ind,nw_ang,fileout)
      
      !get NTO 
      call job_get_nto(ntot,nocc,nvir,nlin,nocc_sym,nvir_sym,movec,fileout)
      
      !deallocate memory
      if (allocated(movec)) deallocate(movec)
      if (allocated(nw_nuc)) deallocate(nw_nuc)
      if (allocated(nw_ind)) deallocate(nw_ind)
      if (allocated(nw_ang)) deallocate(nw_ang)
      call info_mem_deallocation
      call movec_mem_deallocation
      call job_mem_deallocation

   elseif (bastype==2) then
      call read_aces2_info(1,ind1,yesno,fileout,'NBASCART')
      ncart=ind1(1)
      call read_aces2_info(1,ind1,yesno,fileout,'NBASTOT')
      nsph=ind1(1)
      write(fileout,'(a,i5)') 'NBASCART: ',ncart
      write(fileout,'(a,i5)') 'NBASTOT:  ',nsph
   
      natom = 10
      mxang = 4

      !make spherical(nwchem)->cartesian(aces2) matrix: trmat
      allocate(trmat(nsph,ncart))
      call job_nwchem_mo_transform(natom,mxang,ncart,nsph,trmat,fileout)

      return !FIXME

      !get mo
      allocate(movec(nsph,nsph))
      call movec_read_mo(nsph,movec,fileout,uhf)

      !new movec
      call aovec_rotate_sph(nsph,ncart,movec,trmat,fileout)

      if (allocated(movec)) deallocate(movec)
      if (allocated(trmat)) deallocate(trmat)
   endif
end subroutine job_nwchem_nto


subroutine job_rotate_ao(ntot,nocc,nvir,movec,aolab, &
              orbsym_occ,orbsym_vir,nw_nuc,nw_ind,nw_ang,fileout)
   use mod_print
   use mod_movec
   use mod_civec
   use mod_get_aces2
   implicit none
   integer ntot,nocc,nvir,fileout
   double precision movec(ntot,ntot)
   double precision MatTr1(ntot,ntot)
   double precision MatTr2(ntot,ntot)

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
   integer ao_ind(ntot,3),ao_ang(ntot,3)

   logical lnosym

!  lnosym = .false.
   lnosym = .true.

!  get mo (from nwchem), overlap (from aces2) and calculate C'SC
   allocate(MatS(ntot,ntot))
!  call get_aces2_overlap(ntot,MatS,fileout,1)

!  ao_nuc/ao_str/ao_ind/ao_ang
   write(fileout,'(/,a)') 'Print read BASINFO'
   call get_aces2_basinfo(ntot,ao_nuc,ao_str,ao_ind,ao_ang, &
        'ACES2_BASINFO',fileout)

!  read S, AOSO, SOAO matrix from ACES2
   call get_aces2_matrix1(ntot,MatS,'alice_aces2_overlap.txt',fileout)
   write(fileout,'(/,a)') 'Print Overlap (original)'
   call print_r2mat(MatS,ntot,ntot,5,fileout,1)

!  call get_aces2_matrix1(ntot,MatTr1,'alice_aces2_mat_aoso.txt',fileout)
!  write(fileout,'(/,a)') 'Print transformation matrix (AOSO)'
!  call print_r2mat(MatTr1,ntot,ntot,5,fileout,2)

   call get_aces2_matrix1(ntot,MatTr2,'alice_aces2_mat_soao.txt',fileout)
   write(fileout,'(/,a)') 'Print transformation matrix (SOAO)'
   call print_r2mat(MatTr2,ntot,ntot,5,fileout,2)

!  Rotato AO ordering (aorot(NWCHEM/ACES2) -> movec, aolab)
   call movec_match_aoind(ntot,ao_ind,ao_ang,ao_str,fileout,1)
!  aorot(new)
   call movec_aorot_nw2ac(ntot,nw_ind,nw_ang,nw_nuc, &
        ao_ind,ao_ang,ao_nuc,aorot,fileout)
!  call movec_aorot_nw2ac2(ntot,nw_ind,nw_ang,nw_nuc, &
!       ao_ind,ao_ang,ao_nuc,aorot,fileout)

!  aorot(old)
!  call movec_get_aorot(ntot,nw_nuc,nw_ind,nw_ang,aorot, &
!          ao_nuc,ao_str,ao_ind,ao_ang,fileout)

!  rotate movec/aolab/nuc_name
   call movec_ao_rotate1(ntot,movec,aorot,aolab,nw_nuc,fileout,0)
   if (lnosym) then
      ! get ao_norm 
      call movec_ao_norm1(ntot,nw_ang,ao_norm,fileout)
!     movec -> mavec(aorot)
      call movec_ao_norm(ntot,movec,MatS,aolab,aorot,ao_norm,fileout,2)
   endif

   write(fileout,'(/,a)') "* Print AO rotated MO vector1 (NWChem w/AO lable)"
   call print_mo(movec,ntot,ntot,ntot,nocc,nvir, &
        orbsym_occ,orbsym_vir,nw_nuc,aolab,5,fileout)

   ! Rotate AO ordering (SOAO)
   call movec_ao_rotSOAO(ntot,movec,MatTr2,aolab,nw_nuc,fileout)

   ! test
   !call movec_rotate_soao(ntot,movec,MatTr1,MatTr2,fileout)
   !call movec_rotate_soao(ntot,MatS,MatTr1,MatTr2,fileout)

   write(fileout,'(/,a)') "* Print AO rotated MO vector2 (NWChem w/AO lable)"
   call print_mo(movec,ntot,ntot,ntot,nocc,nvir, &
        orbsym_occ,orbsym_vir,nw_nuc,aolab,5,fileout)

   call print_overlap(ntot,movec,movec,MatS,fileout,'CSC')
end subroutine job_rotate_ao 


subroutine job_get_nto(ntot,nocc,nvir,nlin,nocc_sym,nvir_sym,movec,fileout)
   use mod_print
   use mod_info
   use mod_movec
   use mod_civec
   use mod_get_aces2
   implicit none
   integer ntot,nocc,nvir,nlin,fileout
   integer nocc_sym(8),nvir_sym(8)
   double precision movec(ntot,ntot)
   double precision, allocatable :: mo_nto(:,:) 
   character(5) molab(8)
   integer nirrep,occupya(8),numbasir(8),nlindep(8)

   !get ACES2 info 
   call get_aces2_info(nirrep,occupya,numbasir,nlindep,molab,fileout)

   call civec_read_tmat(ntot,nocc,nvir,nlin,orbeng_occ,orbeng_vir,fileout)
   call civec_tmat_rotate(nocc,nvir,rot_occ,rot_vir,fileout)
   !test
!  call civec_get_transvec2(tmat2,nocc,nvir,nocc_sym,nvir_sym, &
!       orbind2_occ,orbind2_vir,orbsym2_occ,orbsym2_vir,fileout)

   call civec_get_transvec(nocc,nvir,nocc_sym,nvir_sym, &
        orbind2_occ,orbind2_vir,orbsym2_occ,orbsym2_vir,fileout) 

   !get NTO and write into file
   allocate(mo_nto(ntot,ntot))
   call movec_get_nto(ntot,nocc,nvir,ntoval1,ntoval2,ntovec1,ntovec2, &
        movec,mo_nto,fileout)
   call movec_write_mo_aces2(ntot,nocc,nvir,nocc_sym,nvir_sym,mo_nto, &
          nirrep,occupya,numbasir,fileout)
   deallocate(mo_nto)
end subroutine job_get_nto 


subroutine job_mem_deallocation
   implicit none
   if(allocated(MatS)) deallocate(MatS)
end subroutine job_mem_deallocation

end module mod_job
