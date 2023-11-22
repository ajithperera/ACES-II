module mod_info
   integer, allocatable :: orbind_occ(:)
   integer, allocatable :: orbind_vir(:)
   integer, allocatable :: orbind2_occ(:)
   integer, allocatable :: orbind2_vir(:)
   integer, allocatable :: rot_occ(:)
   integer, allocatable :: rot_vir(:)

   character(len=5), allocatable :: orbsym_occ(:)
   character(len=5), allocatable :: orbsym_vir(:)
   character(len=5), allocatable :: orbsym2_occ(:)
   character(len=5), allocatable :: orbsym2_vir(:)

   double precision, allocatable :: orbeng_occ(:)
   double precision, allocatable :: orbeng_vir(:)

   character(len=10),allocatable :: aolab(:)
contains

subroutine info_read_mo(ntot,nocc,nvir,nlin,nirrep,nlindep,fileout,lprint)
   implicit none
   integer fileinp,fileout,nirrep,nlindep(8)
   integer ntot,nocc,nvir,nlin,ind,i,j,a,itemp
   logical lprint
   character(80) string
   character(5)  nw_symlab(8)
   integer nw_symnum(8),isymnum(8),nlin_sym

   fileinp = 21
   open(unit=fileinp,file='nwchem.info',status='old',position='rewind')
   read(fileinp,*) string, ntot
   read(fileinp,*) string, nocc
   read(fileinp,*) string, nvir
   read(fileinp,*) string, nlin  !number of linear dependent orbitals
   read(fileinp,*) string, nirrep 
   read(fileinp,*) string, (nw_symlab(i),i=1,nirrep)
   read(fileinp,*) string, (nw_symnum(i),i=1,nirrep)
   read(fileinp,*) string

   allocate(orbind_occ(nocc))
   allocate(orbsym_occ(nocc))
   allocate(orbeng_occ(nocc))
   allocate(orbind_vir(nvir))
   allocate(orbsym_vir(nvir))
   allocate(orbeng_vir(nvir))

   isymnum(1:8)=0
   do i = 1, nocc
      read(fileinp,*) orbind_occ(i), orbsym_occ(i), itemp, orbeng_occ(i)
      do j = 1, nirrep
         if(orbsym_occ(i)==nw_symlab(j)) isymnum(j)=isymnum(j)+1
      enddo
   enddo
   do i = 1, nvir-nlin
      read(fileinp,*) orbind_vir(i), orbsym_vir(i), itemp, orbeng_vir(i)
      do j = 1, nirrep
         if(orbsym_vir(i)==nw_symlab(j)) isymnum(j)=isymnum(j)+1
      enddo
   enddo

   close(fileinp)

   ind = nvir-nlin
   do i = 1, nirrep
      nlin_sym = nw_symnum(i) - isymnum(i)
      do j = 1, nlin_sym
         ind = ind + 1 
         orbind_vir(ind) = nocc+ind
         orbsym_vir(ind) = nw_symlab(i)
         orbeng_vir(ind) = 10000.0D0 
      enddo
   enddo

   if (lprint) then
      write(fileout,'(/,a)') "*info_read_mo: Print symmetry lable (occupied)"
      write(fileout,'(5(i5,a,a5))') (orbind_occ(i),":",orbsym_occ(i),i=1,nocc)
      write(fileout,'(/,a)') "*info_read_mo: Print symmetry lable (virtual)"
      !write(fileout,'(5(i5,a,a5))') (orbind_vir(i),":",orbsym_vir(i),i=1,nvir-nlin)
      write(fileout,'(5(i5,a,a5))') (orbind_vir(i),":",orbsym_vir(i),i=1,nvir)
      write(fileout,*)
      write(fileout,'(a,i3)') "* Number of linear dependent orbitals = ",nlin
      write(fileout,'(a,8a5)') "* Symmetry labels  :    ",(nw_symlab(i),i=1,nirrep)
      write(fileout,'(a,8i5)') "* Orbitals for sym.: ",(nw_symnum(i),i=1,nirrep)
      write(fileout,'(a,8i5)') "* Orb. in final MO : ",(isymnum(i),i=1,nirrep)
   endif
end subroutine info_read_mo



subroutine info_read_acesgeo(natom,geoind,fileout)
   implicit none
   integer natom,geoind(natom,3)
   integer fileinp,fileout
   integer i,j,ind,num(3)
   logical lread
   character(2) str(2)
   character(80) string
   double precision geo(3)

   fileinp = 21
   lread=.false.
   ind=0
   open(unit=fileinp,file='ACES2_BASINFO',status='old',position='rewind')
   do while(.true.)
      read(fileinp,'(a80)',end=10) string

      if (string(1:1)=='*') then
         if (string(1:10)=='* geometry') then
            lread=.true.
         else
            lread=.false.
         endif 
      endif

      if (lread .and. string(1:1)/='*') then
         read(string,*) str(1),num(1),str(2),num(2),num(3),geo(1),geo(2),geo(3)
         ind=ind+1
         geoind(ind,1)=num(1)
         geoind(ind,2)=num(2)
         geoind(ind,3)=num(3)
      endif
   enddo

10 continue

   write(fileout,'(/,a)') '* geometry'
   do i = 1, natom
      write(fileout,'(3i5)') (geoind(i,j),j=1,3)
   enddo

   close(fileinp)
end subroutine info_read_acesgeo


subroutine info_read_aolab(ntot,nw_nuc,nw_ind,nw_ang,fileout)
!for cartesian basis functions
   implicit none
   integer ntot,fileout
   integer fileinp,i,j,nnuc,maxind,itemp
   character(len=10) str_ang_read
   character(len=1) str_ang(7)
   character(len=3) nw_nuc(ntot)
   integer nw_ind(ntot,3)
   integer nw_ang(ntot,3)
   data (str_ang(i),i=1,7) /'s','p','d','f','g','h','i'/

   allocate(aolab(ntot))
   
   fileinp = 21
   open(unit=fileinp,file='nwchem.aolab',status='old',position='rewind')
   nnuc = 0
   do i = 1, ntot
      read(fileinp,*) j,nw_ind(i,1),nw_nuc(i),str_ang_read
      aolab(i) = str_ang_read
      if (nw_ind(i,1)>nnuc) nnuc=nw_ind(i,1) 

      ! nw_ang(:,1) = s~i -> 1~7
      do j = 1, 7
         if (str_ang_read(1:1)==str_ang(j)) then
           nw_ang(i,1) = j
         endif
      enddo

      ! nw_ang(:,2) = 100x+10y+z
      nw_ang(i,2) = 0
      do j = 2, 10
         if (str_ang_read(j:j)=='x') then
            nw_ang(i,2) = nw_ang(i,2) + 100
         elseif (str_ang_read(j:j)=='y') then
            nw_ang(i,2) = nw_ang(i,2) + 10
         elseif (str_ang_read(j:j)=='z') then
            nw_ang(i,2) = nw_ang(i,2) + 1
         endif
      enddo
   enddo

   ! nw_ang(:,3)
   nw_ang(1,3) = 1
   do i = 2, ntot
      nw_ang(i,3) = 0
      maxind = 0
      do j = 1, i-1
         if (nw_ind(i,1)==nw_ind(j,1)) then
            if (nw_ang(i,1)==nw_ang(j,1)) then         
            if (nw_ang(i,2)==nw_ang(j,2)) then
               if (nw_ang(j,3)>maxind) maxind=nw_ang(j,3)
            endif
            endif
         endif
      enddo
      nw_ang(i,3)=maxind+1 
   enddo

   do i = 1, ntot
      write(fileout,'(a,3i5)') 'nw_ang',nw_ang(i,1),nw_ang(i,2),nw_ang(i,3)
   enddo

   close(fileinp)
end subroutine info_read_aolab


subroutine info_symblk(nocc,nvir,nlin,nocc_sym,nvir_sym,nirrep,molab,fileout)
   implicit none
   integer      nocc,nvir,nlin,fileout
   character(5) symlab(8),molab(8)
   integer      i,j,isym,iocc,ivir,ind,ind_occ(nocc),ind_vir(nvir)
   integer      nocc_sym(8),nvir_sym(8)
   integer  nirrep
   do i = 1, 8
      symlab(i) = "     "
      nocc_sym(i) = 0
      nvir_sym(i) = 0
   enddo

   do i = 1, nirrep
      symlab(i)=lower(molab(i))
   enddo
   write(fileout,'(/,a)') 'MO label ordering (ACES2)' 
   write(fileout,'(8(i2,a3,a5))') (i,' : ',symlab(i),i=1,nirrep) 

   do i = 1, nocc
      do j = 1, 8
         if (symlab(j)=="     ") then
            symlab(j)=orbsym_occ(i)
            ind_occ(i)=j
            nocc_sym(j) = nocc_sym(j) + 1
            exit
         elseif (symlab(j)==orbsym_occ(i)) then
            ind_occ(i)=j
            nocc_sym(j) = nocc_sym(j) + 1
            exit
         endif
      enddo
   enddo

   do i = 1, nvir
      do j = 1, 8
         if (symlab(j)=="     ") then
            symlab(j)=orbsym_vir(i)
            ind_vir(i)=j
            nvir_sym(j) = nvir_sym(j) + 1
            exit
         elseif (symlab(j)==orbsym_vir(i)) then
            ind_vir(i)=j
            nvir_sym(j) = nvir_sym(j) + 1
            exit
         endif
      enddo
   enddo

!  rotate index (rot_occ) for occupied
   allocate(rot_occ(nocc))
   ind = 0
   do isym = 1, 8
      do iocc = 1, nocc
         if (ind_occ(iocc)==isym) then
            ind = ind + 1
            rot_occ(ind) = iocc
         endif
      enddo
   enddo

!  rotate index (rot_vir) for virtual
   allocate(rot_vir(nvir))
   ind = 0
   do isym = 1, 8
      do ivir = 1, nvir
         if (ind_vir(ivir)==isym) then
            ind = ind + 1
            rot_vir(ind) = ivir
         endif
      enddo
   enddo

!  print
   write(fileout,'(/,a)') "*info_symblk: Print symmetry label (occupied)"
   write(fileout,'(5(i5,a,a5))') (ind_occ(i),":",orbsym_occ(i),i=1,nocc)

   write(fileout,'(/,a)') "*info_symblk: Print symmetry label (virtual)"
   !write(fileout,'(5(i5,a,a5))') (ind_vir(i),":",orbsym_vir(i),i=1,nvir-nlin)
   write(fileout,'(5(i5,a,a5))') (ind_vir(i),":",orbsym_vir(i),i=1,nvir)
   write(fileout,'(a,i3)') "  Number of linear dependent orbitals = ",nlin

   write(fileout,'(/,a)') "*info_symblk: Print symmetry ordering (occupied)"
   write(fileout,'(10(i5))') (rot_occ(i),i=1,nocc)

   write(fileout,'(/,a)') "*info_symblk: Print symmetry ordering (virtual)"
   !write(fileout,'(10(i5))') (rot_vir(i),i=1,nvir-nlin)
   write(fileout,'(10(i5))') (rot_vir(i),i=1,nvir)
   write(fileout,'(a,i3)') "  Number of linear dependent orbitals = ",nlin
end subroutine info_symblk


subroutine info_get_rotindex(nocc,nvir)
   implicit none
   integer nocc,nvir
   integer i,a

   allocate(orbind2_occ(nocc))
   allocate(orbsym2_occ(nocc))
   do i = 1, nocc
      orbind2_occ(i) = orbind_occ(rot_occ(i))
      orbsym2_occ(i) = orbsym_occ(rot_occ(i))
   enddo

   allocate(orbind2_vir(nvir))
   allocate(orbsym2_vir(nvir))
   do i = 1, nvir
      orbind2_vir(i) = orbind_vir(rot_vir(i))
      orbsym2_vir(i) = orbsym_vir(rot_vir(i))
   enddo
end subroutine info_get_rotindex

subroutine info_mem_deallocation
   implicit none
   if (allocated(orbind_occ)) deallocate(orbind_occ)
   if (allocated(orbind_vir)) deallocate(orbind_vir)
   if (allocated(orbind2_occ)) deallocate(orbind2_occ)
   if (allocated(orbind2_vir)) deallocate(orbind2_vir)
   if (allocated(orbsym_occ)) deallocate(orbsym_occ)
   if (allocated(orbsym_vir)) deallocate(orbsym_vir)
   if (allocated(orbsym2_occ)) deallocate(orbsym2_occ)
   if (allocated(orbsym2_vir)) deallocate(orbsym2_vir)
   if (allocated(orbeng_occ)) deallocate(orbeng_occ)
   if (allocated(orbeng_vir)) deallocate(orbeng_vir)
   if (allocated(rot_occ)) deallocate(rot_occ)
   if (allocated(rot_vir)) deallocate(rot_vir)
   if (allocated(aolab)) deallocate(aolab)
end subroutine info_mem_deallocation

FUNCTION Lower(s1)  RESULT (s2)
   CHARACTER(*)       :: s1
   CHARACTER(LEN(s1)) :: s2
   CHARACTER          :: ch
   INTEGER,PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
   INTEGER            :: i

   DO i = 1,LEN(s1)
      ch = s1(i:i)
      IF (ch >= 'A'.AND.ch <= 'Z') ch = CHAR(ICHAR(ch)-DUC)
      s2(i:i) = ch
   END DO
END FUNCTION Lower

end module mod_info
