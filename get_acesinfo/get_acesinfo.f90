program get_acesinfo
   implicit none
   integer fileinp,fileout
   character(len=80) keyword
   integer iuhf
!  COMMON BLOCKS

integer :: icore(1)
common / / icore

! istart.com : begin
      integer         i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
! istart.com : end
! flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
! flags.com : end
! flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
! flags2.com : end

   fileinp = 11
   fileout = 6 
!  initialize
   call aces_init(icore,i0,icrsiz, iuhf, .true.)

!  read input file (a2info.inp)
   open (unit=fileinp,file='a2info.inp',status='old') 
   open (unit=fileout,file='a2info.out',status='unknown') 

   write (fileout,'(a,/)') '@ get_acesinfo : Start program'

   do while(.true.)
      read (fileinp,'(a)',end=90) keyword
      !write (fileout,'(x,a)') keyword
      call read_JOBARC(trim(keyword),fileout)

      call get_moinfo(trim(keyword),fileout)
   end do

90 continue
   write (fileout,'(/,a)') '@ get_acesinfo : End of program'

   close (fileinp)
   close (fileout)

   call aces_fin
end program get_acesinfo
