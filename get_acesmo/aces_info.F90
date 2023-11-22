subroutine get_aces2_matrix1(ntot,Mat1,fileread,fileout)
   implicit none
   integer ntot,fileout
   character(len=*) fileread
   double precision Mat1(ntot,ntot)
   integer fileinp,i,j,ind,itmp,jtmp

   fileinp = 21
   open(unit=fileinp,file=fileread,status='old',position='rewind')
   do j = 1, ntot
      do i = 1, ntot
         read (fileinp,*) jtmp,itmp,Mat1(i,j)
      enddo
   enddo
   close(fileinp)
end subroutine get_aces2_matrix1


subroutine get_aces2_info(nirrep,occupya,occupyb,numbasir,nlindep,molab,fileout)
   implicit none
   character(80) string
   integer nirrep,occupya(8),occupyb(8),numbasir(8),nlindep(8)
   character(5) molab(8)
   integer i,fileinp,fileout

   fileinp = 21
   open(unit=fileinp,file="a2info.dat",status='old',position='rewind')
   do while (.true.)
      read(fileinp,'(A80)',end=90) string
      if (string(1:7)=='*NIRREP') then 
         read(string(11:80),*) nirrep
         write(fileout,'(a,i5)') 'ACES2,NIRREP : ',nirrep
      elseif (string(1:8)=='*OCCUPYA') then 
         read(string(11:80),*) (occupya(i),i=1,nirrep)
         write(fileout,'(a,8i5)') 'ACES2,OCCUPYA : ',(occupya(i),i=1,nirrep)
      elseif (string(1:8)=='*OCCUPYB') then 
         read(string(11:80),*) (occupyb(i),i=1,nirrep)
         write(fileout,'(a,8i5)') 'ACES2,OCCUPYB : ',(occupyb(i),i=1,nirrep)
      elseif (string(1:9)=='*NUMBASIR') then 
         read(string(11:80),*) (numbasir(i),i=1,nirrep)
         write(fileout,'(a,8i5)') 'ACES2,NUMBASIR : ',(numbasir(i),i=1,nirrep)
      elseif (string(1:6)=='*MOLAB') then 
         if (nirrep==1) then
            molab(1)='A1'
         else
            read(string(11:80),*) (molab(i),i=1,nirrep)
         endif
      elseif (string(1:8)=='*NLINDEP') then 
         read(string(11:80),*) (nlindep(i),i=1,nirrep)

      endif
   enddo
90 continue
   close(fileinp)
   return
end subroutine get_aces2_info


subroutine get_aces2_basinfo(ntot,ao_nuc,ao_str,ao_ind,ao_ang, &
           fileread,fileout)
   implicit none
   character(len=*) fileread
   character(80) string
   integer ntot
   character(4), allocatable :: nuc(:)
   character(1), allocatable :: nuctype(:)
   integer, allocatable :: nucnum(:,:)
   double precision, allocatable :: geo(:,:)
   integer num,i,j,k
   integer ind1,ind2,ind3,ind4,ind5,ind6,ind7
   character(2) str1
   character(5) str2
   integer fileinp,fileout
   character(2) ao_nuc(ntot) 
   character(5) ao_str(ntot) 
   integer ao_ind(ntot,3),ao_ang(ntot,3)

   fileinp = 21
   open(unit=fileinp,file=fileread,status='old',position='rewind')
   num = 0
   read(fileinp,'(A80)',end=90) string
   do while (.true.)
      read(fileinp,'(A80)',end=90) string
      if (string(1:1)=='*') exit
      num = num + 1
   enddo

   !geometry
   write(fileout,'(a,x,i3)') 'geometry lines',num
   allocate(nuc(num))
   allocate(nuctype(num))
   allocate(nucnum(num,3))
   allocate(geo(num,3))

   rewind(fileinp)
   read(fileinp,'(a80)',end=90) string

   i = 0
   do while (.true.) 
      read(fileinp,'(a80)',end=90) string
      if (string(1:1)=='*') exit      
      i = i + 1
      read(string,*) &
          nuctype(i),nucnum(i,1),nuc(i),nucnum(i,2),nucnum(i,3), &
          geo(i,1),geo(i,2),geo(i,3) 
   enddo
   Print*, nuctype(i),nucnum(i,1),nuc(i),nucnum(i,2),nucnum(i,3), &
          geo(i,1),geo(i,2),geo(i,3)

   !AO list
   !write(fileout,'(a)') 'reading AO list'
   i = 0
   do while (.true.) 
      read(fileinp,'(a80)',end=90) string
      if (string(1:1)=='*') exit      
      read(string,*) str1,ind1,ind2,ind3,ind4,str2,ind5,ind6,ind7 

      i = i + 1
      ao_nuc(i) = str1
      ao_str(i) = str2
      !ao_ind(:,1)
      do k = 1, num
         if (str1==nuc(k) .and. ind1==nucnum(k,2) .and. ind2==nucnum(k,3)) then
            ao_ind(i,1) = nucnum(k,1) 
         endif
      enddo
      write(fileout,'(a,2i5)')'ao_ind1',i,ao_ind(i,1)
      ao_ind(i,2) = ind1  !atom number(sym-indep)
      ao_ind(i,3) = ind2  !atom number(all)
      ao_ang(i,1) = ind3  !ang.mom(1,3,6,10,...)
      ao_ang(i,2) = ind4  !ang.mom order
                  ! ind5  !No. primitive
                  ! ind6  !No. contraction
      ao_ang(i,3) = ind7  !contraction order
   enddo
90 continue
   close(fileinp)

!  AO list
!  ao_nuc(:) : name of atoms
!  ao_ind(:,1~3) : order of atom / sym. indep.-atom / parity?  
!  ao_ang(:,1~2) : total possible l,m,n combination (x^l y^m z^n) / ordering
!  ao_str(:)
   write(fileout,'(/,a)') 'Print AO list'
   do i = 1, ntot
      write(fileout,'(i3,x,a2,3i3,x,a,x,2i3,x,a5)') i,ao_nuc(i),(ao_ind(i,j),j=1,3), &
           ":",(ao_ang(i,j),j=1,2),ao_str(i)
   enddo

   deallocate(nuc)
   deallocate(nuctype)
   deallocate(nucnum)
   deallocate(geo)
end subroutine get_aces2_basinfo

