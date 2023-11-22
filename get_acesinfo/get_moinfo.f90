subroutine get_moinfo(string,fileout)
   implicit none
   logical yesno
   character(len=*) string
   character(8), allocatable :: molab(:)
   character(8)  lab(8)
   integer,  allocatable :: ivmlsym(:)
   integer fileread,filework,fileout
   integer nirrep,occupya(8),numbasir(8),nbfirr(8),nlindep(8)
   integer i,j,nbas,iref

   fileread = 21
   filework = 25

   if (string.eq.'ACES2.MOINFO') then

      inquire(file='a2info.dat',exist=yesno)
      if (yesno) then
         write (fileout,'(a)') ' a2info.dat exists and will be deleted.'
         open(filework,file='a2info.dat',status='old',form='formatted')
         close(unit=filework,status='delete')
      endif
      open(filework,file='a2info.dat',status='new',form='formatted')

      !NIRREP
      call getrec(fileread,'JOBARC','NIRREP',1,nirrep)
      write (filework,'(a,i5)') "*NIRREP  :",nirrep
      
      !OCCUPYA : Number of occupancy per irrep.
      call getrec(fileread,'JOBARC','OCCUPYA',nirrep,occupya)
      write (filework,'(a,8i5)') "*OCCUPYA :",(occupya(i),i=1,nirrep)

      call getrec(fileread,'JOBARC','UHFRHF',1,iref)
!     write(fileout,'(a,i1)') 'UHFRHF = ',iref
      if (iref==1) then
         call getrec(fileread,'JOBARC','OCCUPYB',nirrep,occupya)
         write (filework,'(a,8i5)') "*OCCUPYB :",(occupya(i),i=1,nirrep)
      endif

      !NUMBASIR : Number of basis functions per irrep.
      call getrec(fileread,'JOBARC','NUMBASIR',nirrep,numbasir)
      write (filework,'(a,8i5)') "*NUMBASIR:",(numbasir(i),i=1,nirrep)

      !NBFIRR (=NUMBASIR ?)
!     call getrec(fileread,'JOBARC','NBFIRR',nirrep,nbfirr)
      call getrec(fileread,'JOBARC','NUMBASIR',nirrep,nbfirr)
      write (filework,'(a,8i5)') "*NUMBASIR:",(nbfirr(i),i=1,nirrep)

      !NLINDEP : Number of linear dependent orbitals per irrep.
      call getrec(fileread,'JOBARC','NLINDEP',8,nlindep)
      write (filework,'(a,8i5)') "*NLINDEP :",(nlindep(i),i=1,nirrep)

      call getrec(fileread,'JOBARC','NBASCART',1,nbas)  !Num0=NAO
      write (filework,'(a,i5)') "*NBASCART:",nbas

      call getrec(fileread,'JOBARC','NBASTOT',1,nbas)  !Num0=NAO
      write (filework,'(a,i5)') "*NBASTOT :",nbas

      !GET MO ordering
      allocate (molab(nbas))
      allocate (ivmlsym(nbas))
      call GETREC(fileread,'JOBARC','EVCSYMAC',nbas,molab)
      call GETREC(fileread,'JOBARC','IVMLSYM',nbas,ivmlsym)
      
      do i = 1, nirrep
         yesno = .false.
         do j = 1, nbas
            if (ivmlsym(j)==i .and. yesno==.false.) then
               lab(i)=molab(j)  
               yesno = .true.
            endif
         enddo
      enddo
      write(filework,'(a,8(x,a4))') "*MOLAB   :",(lab(i)(1:4),i=1,nirrep)
      deallocate (molab)
      deallocate (ivmlsym)

      close(filework)
   endif

end subroutine get_moinfo
