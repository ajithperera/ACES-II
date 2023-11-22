!subroutine read_JOBARC(string,Num0,Num1,Val0,Val1,fileout)
subroutine read_JOBARC(string,fileout)
   implicit none
   character(len=*) string
   character(len=80) str0
   integer,          allocatable :: Num1(:)
   character(8),     allocatable :: Str8(:)
   double precision, allocatable :: Val1(:)
   integer fileread, fileout
   integer IINTFP,Num0,Dim0,Dim1,Dim2
   double precision Val0
   integer i,j,k
!  integer              :: Num0
!  integer, allocatable :: Num1(:)
!  double precision              :: Val0
!  double precision, allocatable :: Val1(:)

   IINTFP = 1  !64bit=1, 32bit=2
   fileread = 21

   if (string.eq.'NATOMS') then
      call getrec(fileread,'JOBARC','NATOMS',1,Num0)
      write (fileout,'(a,i5)') '* NATOMS : Number of atoms = ',Num0

   else if (string.eq.'FULLAOSO') then
!     12.
      call getrec(fileread,'JOBARC','NBASTOT',1,Num0)  !NAO
      Dim0 = Num0*Num0   !NAO*NAO
      allocate (Val1(Dim0))
      call getrec(fileread,'JOBARC','FULLAOSO',Dim0,Val1)
      write (fileout,'(a,i5)') "* FULLAOSO : Transformation matrix for SO = ",Num0
!     call print_rmat(Val1,Dim0,Num0,Num0,fileout,2)
      call print_r2mat(Val1,Num0,Num0,5,fileout,2)
      call write_aces_r2array(Val1,Dim0,Num0,Num0, &
           'a2info_mat_aoso.dat',fileout)

   else if (string.eq.'FULLSOAO') then
!     13. 
      call getrec(fileread,'JOBARC','NBASTOT',1,Num0)  !NAO
      Dim0 = Num0*Num0   !NAO*NAO
      allocate (Val1(Dim0))
      call getrec(fileread,'JOBARC','FULLSOAO',Dim0,Val1)
      write (fileout,'(a,i5)') "* FULLSOAO : Transformation matrix for SO = ",Num0
!     call print_rmat(Val1,Dim0,Num0,Num0,fileout,2)
      call print_r2mat(Val1,Num0,Num0,5,fileout,2)
      call write_aces_r2array(Val1,Dim0,Num0,Num0, & 
           'a2info_mat_soao.dat',fileout)

   else if (string.eq.'NBASTOT') then
!     27. Number of symmetry-adapted orbitals in the basis
      call getrec(fileread,'JOBARC','NBASTOT',1,Num0)  !NAO? NSO?
      write (fileout,'(a,i5)') "* NBASTOT : Number of atomic orbitals = ",Num0

   else if (string.eq.'NOCCORB') then
!     29. The number of occupied molecular orbitals
      allocate (Num1(2))
      call getrec(fileread,'JOBARC','NOCCORB',2,Num1)  !NOCC
      write (fileout,'(a,2i5)') "* NOCCORB : Number of occupied orbitals = ", &
            Num1(1),Num1(2)

   else if (string.eq.'NUMBASIR') then
!     31. Number of basis functions per irrep.      
      call getrec(fileread,'JOBARC','NIRREP',1,Num0)
      allocate (Num1(8))
      Num1(1:8) = 0
      call getrec(fileread,'JOBARC','NUMBASIR',Num0,Num1)
      write (fileout,'(a,8i5)') "* NUMBASIR : Number of basis functions per irrep.= ", &
            (Num1(I),I=1,Num0)

   else if (string.eq.'NVRTORB') then
!     33. The number of virtual molecular orbitals       
      allocate (Num1(2))
      call getrec(fileread,'JOBARC','NVRTORB',2,Num1)  !NVRT
      write (fileout,'(a,2i5)') "* NVRTORB : Number of virtual orbitals = ",Num1(1),Num1(2)

   else if (string.eq.'OCCUPYA') then
!     35. The number of occupancy per irrep.
      allocate (Num1(8))
      Num1(1:8) = 0
      call getrec(fileread,'JOBARC','NIRREP',1,Num0) 
      call getrec(fileread,'JOBARC','OCCUPYA',Num0,Num1)  !NVRT
      write (fileout,'(a,8i5)') "* OCCUPYA : Number of occupancy per irrep. = ", &
            (Num1(i),i=1,Num0)

   else if (string.eq.'SCFENEG') then
      !call getrec(fileread,'JOBARC','SCFENEG',IINTFP,Val0)  !Val0=ESCF

   else if (string.eq.'NUCREP') then
      call getrec(fileread,'JOBARC','NUCREP',1,Val0)
      write (fileout,'(a,f15.7)') '* NUCREP : Nuclear repulsion energy = ',Val0

   else if (string.eq.'SCFEVCA0') then
!     45. MO for alpha spin (SCF orbital numbering)
      allocate (Num1(6))
      call getrec(fileread,'JOBARC','NBASTOT',1,Num1(1))  !NAO
      call getrec(fileread,'JOBARC','NOCCORB',2,Num1(2))  !NOCC
      call getrec(fileread,'JOBARC','NVRTORB',2,Num1(4))  !NVRT
      Num1(6) = Num1(2) + Num1(4)  !Num1(6)=NOCC+NVRT
      !write (fileout,'(6i5)') (Num1(i),i=1,6)
      Dim0 = Num1(1)*Num1(6)   !NAO*NMO
      allocate (Val1(Dim0))
      write (fileout,'(2a)') "* SCFEVCA0 : MO for alpha spin",  &
                             " (SCF orbital numbering)"
      call getrec(fileread,'JOBARC','SCFEVCA0',Dim0,Val1)
      call print_rmat(Val1,dim0,Num1(1),Num1(6),fileout,2)


   else if (string.eq.'SCFEVCAS') then
!     MO for alpha spin (SCF orbital numbering)
      allocate (Num1(6))
      call getrec(fileread,'JOBARC','NBASTOT',1,Num1(1))  !NAO
      call getrec(fileread,'JOBARC','NOCCORB',2,Num1(2))  !NOCC
      call getrec(fileread,'JOBARC','NVRTORB',2,Num1(4))  !NVRT
      Num1(6) = Num1(2) + Num1(4)  !Num1(6)=NOCC+NVRT
      Dim0 = Num1(1)*Num1(6)   !NAO*NMO
      allocate (Val1(Dim0))
      write (fileout,'(2a)') "* SCFEVCAS : MO for alpha spin",  &
                             " (SCF orbital numbering)"
      call getrec(fileread,'JOBARC','SCFEVCAS',Dim0,Val1)
      call print_rmat(Val1,dim0,Num1(1),Num1(6),fileout,2)


   else if (string.eq.'SCFEVECA') then
!     47. MO for alpha spin  (post-SCF orbital numbering)
      allocate (Num1(6))
      call getrec(fileread,'JOBARC','NBASTOT',1,Num1(1))  !NAO
      call getrec(fileread,'JOBARC','NOCCORB',2,Num1(2))  !NOCC
      call getrec(fileread,'JOBARC','NVRTORB',2,Num1(4))  !NVRT
      Num1(6) = Num1(2) + Num1(4)  !Num1(6)=NOCC+NVRT
      !write (fileout,'(6i5)') (Num1(i),i=1,6)
      Dim0 = Num1(1)*Num1(6)   !NAO*NMO
      allocate (Val1(Dim0))
      write (fileout,'(2a)') "* SCFEVECA : MO for alpha spin",  &
                             " (post-SCF orbital numbering)"
      call getrec(fileread,'JOBARC','SCFEVECA',Dim0,Val1) 
      call print_rmat(Val1,dim0,Num1(1),Num1(6),fileout,2)
      !deallocate (Val1)
      !deallocate (Num1)

   else if (string.eq.'SCFDENSA') then
!     density matrix for ground/excited state/response excited state 
      call getrec(fileread,'JOBARC','NBASTOT',1,Num0)  !Num0=NAO 
      Dim0 = Num0*Num0*IINTFP
      allocate (Val1(Dim0))
      call getrec(fileread,'JOBARC','SCFDENSA',Dim0,Val1)
      write(fileout,'(a)') 'SCFDENSA'
      call print_r2mat(Val1,Num0,Num0,5,fileout,1)
      !deallocate (Val1)

   else if (string.eq.'SCFDENSB') then
!     density matrix for ground/excited state/response excited state 
      call getrec(fileread,'JOBARC','NBASTOT ',1,Num0)  !Num0=NAO 
      Dim0 = Num0*Num0*IINTFP
      allocate (Val1(Dim0))
      call getrec(fileread,'JOBARC','SCFDENSB',Dim0,Val1)
      write(fileout,'(a)') 'SCFDENSB'
      call print_r2mat(Val1,Num0,Num0,5,fileout,1)
      !deallocate (Val1)

   else if (string.eq.'RELDENSA') then
!     density matrix for ground/excited state/response excited state
      call getrec(fileread,'JOBARC','NBASTOT',1,Num0)  !Num0=NAO
      Dim0 = Num0*Num0*IINTFP
      allocate (Val1(Dim0))
      call getrec(fileread,'JOBARC','RELDENSA',Dim0,Val1)
      write (fileout,'(/,a,i7)') "* RELDENSA: Density, alpha ",Num0
      call print_r2mat(Val1,Num0,Num0,5,fileout,1)

   else if (string.eq.'RELDENSB') then
!     density matrix for ground/excited state/response excited state
      call getrec(fileread,'JOBARC','NBASTOT ',1,Num0)  !Num0=NAO
      Dim0 = Num0*Num0*IINTFP
      allocate (Val1(Dim0))
      call getrec(fileread,'JOBARC','RELDENSB',Dim0,Val1)
      write (fileout,'(/,a,i7)') "* RELDENSB: Density, beta ",Num0
      call print_r2mat(Val1,Num0,Num0,5,fileout,1)

   else if (string.eq.'DENSDIFF') then
!     density excited state density different (vee/edens.f)
      call getrec(fileread,'JOBARC','NBASTOT',1,Num0)  !Num0=NAO 
      Dim0 = Num0*Num0*IINTFP
      allocate (Val1(Dim0))
      call getrec(fileread,'JOBARC','DENSDIFF',Dim0,Val1)
      !deallocate (Val1)

   else if (string.eq.'FOCKA') then
!     density excited state density different (vee/edens.f)
      call getrec(fileread,'JOBARC','NBASTOT',1,Num0)  !Num0=NAO 
      Dim0 = Num0*Num0*IINTFP
      allocate (Val1(Dim0))
      call getrec(fileread,'JOBARC','FOCKA',Dim0,Val1)
      write (fileout,'(/,a,i7)') "* FOCKA: Fock matrix alpha ",Num0
      call print_r2mat(Val1,Num0,Num0,5,fileout,1)

   else if (string.eq.'FOCKB') then
!     density excited state density different (vee/edens.f)
      call getrec(fileread,'JOBARC','NBASTOT',1,Num0)  !Num0=NAO 
      Dim0 = Num0*Num0*IINTFP
      allocate (Val1(Dim0))
      call getrec(fileread,'JOBARC','FOCKB',Dim0,Val1)
      write (fileout,'(/,a,i7)') "* FOCKB: Fock matrix beta ",Num0
      call print_r2mat(Val1,Num0,Num0,5,fileout,1)

   else if (string.eq.'TMRIGHT') then
      allocate (Val1(3))
      call getrec(fileread,'JOBARC','TMRIGHT',3*IINTFP,Val1)
      write (fileout,'(a,2i5)') "* TMRIGHT : Transition Moment (right side) ;"
      write (fileout,'(x,3f15.7)') Val1(1),Val1(2),Val1(3)
      !deallocate (Val1)

   else if (string.eq.'AOOVRLAP') then
      call getrec(fileread,'JOBARC','NBASTOT',1,Num0)  !Num0=NAO
      Dim0 = Num0*Num0*IINTFP
      allocate (Val1(Dim0))
      call GETREC(fileread,'JOBARC','AOOVRLAP',Dim0,Val1) 
      write (fileout,'(/,a)') "* Overlap matrix" 
      call print_rmat(Val1,Dim0,Num0,Num0,fileout,2) 

      call write_aces_r2array(Val1,dim0,Num0,Num0, & 
           'a2info_overlap.dat',fileout)

   else if (string.eq.'NBASCART') then
      call getrec(fileread,'JOBARC','NBASCART',1,Num0)  !Num0=NAO
      write (fileout,'(a,i5)') "*NBASCART:",Num0

   else if (string.eq.'CMP2CART') then
      call getrec(fileread,'JOBARC','NBASCART',1,Num0)  !Num0=NAO
      write (fileout,'(a,i5)') "*NBASCART:",Num0
      Dim1=Num0
      call getrec(fileread,'JOBARC','NBASTOT',1,Num0)  !Num0=NAO
      write (fileout,'(a,i5)') "*NBASCART:",Num0
      Dim2=Num0
      Dim0=Dim1*Dim2*IINTFP
      allocate (Val1(Dim0))
      CALL GETREC(20,'JOBARC','CMP2CART',DIM0,Val1)
      write (fileout,'(/,a)') "* CMP2CART matrix"
!     call print_rmat(Val1,Dim0,Dim1,Dim2,fileout,2) 
      call print_r2mat(Val1,Dim1,Dim2,5,fileout,2)

   else if (string.eq.'AO2SO ') then
      call getrec(fileread,'JOBARC','NBASCART',1,Num0)  !Num0=NAO
      write (fileout,'(a,i5)') "*NBASCART:",Num0
      Dim1=Num0
      call getrec(fileread,'JOBARC','NBASTOT',1,Num0)  !Num0=NAO
      write (fileout,'(a,i5)') "*NBASCART:",Num0
      Dim2=Num0
      Dim0=Dim1*Dim2*IINTFP
      allocate (Val1(Dim0))
      CALL GETREC(20,'JOBARC','AO2SO ',DIM0,Val1)
      write (fileout,'(/,a)') "* AO2SO matrix"
!     call print_rmat(Val1,Dim0,Dim1,Dim2,fileout,2)
      call print_r2mat(Val1,Dim1,Dim2,5,fileout,2)
      call write_aces_r2array(Val1,Dim0,Dim1,Dim2, &
           'a2info_mat_ao2so.dat',fileout)

   else if (string.eq.'AO2SOINV') then
      call getrec(fileread,'JOBARC','NBASTOT',1,Num0)  !Num0=NAO
      write (fileout,'(a,i5)') "*NBASCART:",Num0
      Dim1=Num0
      call getrec(fileread,'JOBARC','NBASCART',1,Num0) !Num0=NAO
      write (fileout,'(a,i5)') "*NBASCART:",Num0
      Dim2=Num0
      Dim0=Dim1*Dim2*IINTFP
      allocate (Val1(Dim0))
      CALL GETREC(20,'JOBARC','AO2SOINV',DIM0,Val1)
      write (fileout,'(/,a)') "* AO2SOINV matrix"
!     call print_rmat(Val1,Dim0,Dim1,Dim2,fileout,2)
      call print_r2mat(Val1,Dim1,Dim2,5,fileout,2)
      call write_aces_r2array(Val1,Dim0,Dim1,Dim2, &
           'a2info_mat_ao2so.dat',fileout)

!  Symmetry label, FULLSYM
   else if (string.eq.'EVCSYMAF') then
      call getrec(fileread,'JOBARC','NBASTOT',1,Num0)  !Num0=NAO
      Dim0 = Num0*IINTFP
      allocate (Str8(Dim0))
      call GETREC(fileread,'JOBARC','EVCSYMAF',Dim0,Str8)
      write(fileout,'(/,a)')'* MO label (FULLSYM)'
      write(fileout,*) (Str8(i)(1:4),i=1,Num0)

!  Symmetry label, COMPSYM
   else if (string.eq.'EVCSYMAC') then
      call getrec(fileread,'JOBARC','NBASTOT',1,Num0)  !Num0=NAO
      Dim0 = Num0*IINTFP
      allocate (Str8(Dim0))
      call GETREC(fileread,'JOBARC','EVCSYMAC',Dim0,Str8)
      write(fileout,'(/,a)')'* MO label (COMPSYM)'
      write(fileout,*) (Str8(i)(1:4),i=1,Num0)

      allocate (Num1(Dim0))
      call GETREC(fileread,'JOBARC','IVMLSYM',Dim0,Num1)
      write(fileout,'(/,a)')'* MO label (number)'
      write(fileout,'(1000i4)') (Num1(i),i=1,Num0)

!  Newly implemented keyword (for putmos) 
   else if (string.eq.'NIRREP') then
      call getrec(fileread,'JOBARC','NIRREP',1,Num0) 
      write (fileout,'(a,i5)') "* NIRREP : Number of irrep. = ",Num0 

!  Newly implemented keyword (for putmos) 
   else if (string.eq.'NBFIRR') then
      call getrec(fileread,'JOBARC','NIRREP',1,Num0) 
      allocate (Num1(8))
      Num1(1:8) = 0
      call getrec(fileread,'JOBARC','NBFIRR ',Num0,Num1) 
      write (fileout,'(a,8i5)') "* NBFIRR : Number of basis functions per irrep. = ", &
            (Num1(I),I=1,Num0) 

   else if (string.eq.'DIPOLE_X') then
      call getrec(fileread,'JOBARC','NBASTOT',1,Num0)  !Num0=NAO 
      Dim0 = Num0*Num0*IINTFP
      Dim1 = Num0*(Num0+1)/2
      allocate (Val1(Dim1))
      call getrec(fileread,'JOBARC',string,Dim1,Val1)
      write (fileout,'(/,a,a,a,i7)') "* ",string,": ",Num0
      call print_r2mat_tri(Val1,Num0,Dim1,5,fileout,1)

   end if

   if (allocated(Num1)) deallocate(Num1)
   if (allocated(Val1)) deallocate(Val1)
   if (allocated(Str8)) deallocate(Str8)

end subroutine read_JOBARC

