program alice_nwchem
   use mod_job
   implicit none
   integer fileinp,fileout
   character(len=80) keyword

   fileinp = 7
   fileout = 6

   open (unit=fileinp,file='alice.inp',status='old')
   open (unit=fileout,file='alice.out',status='unknown')

   write (fileout,'(a,/)') '@ alice_nwchem : Start program'
   do while(.true.)
      read (fileinp,'(a)',end=90) keyword
      call job_control(keyword,fileout)
   end do

90 continue
   write (fileout,'(/,a)') '@ alice_nwchem : End of program'

   close (fileinp)
   close (fileout)
end program alice_nwchem
