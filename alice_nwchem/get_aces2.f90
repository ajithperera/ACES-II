module mod_get_aces2
contains

subroutine get_aces2_matrix1(ntot,Mat1,fileread,fileout)
   use mod_print
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


subroutine get_aces2_matrix2(nrow,ncol,Mat1,fileread,fileout)
   use mod_print
   implicit none
   integer nrow,ncol,fileout
   character(len=*) fileread
   double precision Mat1(nrow,ncol)
   integer fileinp,i,j,ind,itmp,jtmp

   fileinp = 21
   open(unit=fileinp,file=fileread,status='old',position='rewind')
   do j = 1, ncol
      do i = 1, nrow
         read (fileinp,*) jtmp,itmp,Mat1(i,j)
      enddo
   enddo
   close(fileinp)
end subroutine get_aces2_matrix2


subroutine get_aces2_overlap(ntot,MatS,fileout,itype)
   use mod_print 
   implicit none
   integer ntot,fileout,itype
   double precision MatS(ntot,ntot)
   double precision MatTr1(ntot,ntot)
   double precision MatTr2(ntot,ntot)

   double precision Mat1(ntot,ntot)

   character(len=80) string
   integer fileinp,i,j,ind,itmp,jtmp,nsym
   double precision temp 

   call get_aces2_matrix1(ntot,MatS,'alice_aces2_overlap.txt',fileout)
   write(fileout,'(/,a)') 'GET_ACES2_OVERLAP: Print Overlap (original)'
   call print_r2mat(MatS,ntot,ntot,5,fileout,2)

   nsym = 1
   if (nsym==1) return

   call get_aces2_matrix1(ntot,MatTr1,'alice_aces2_mat_aoso.txt',fileout)
   write(fileout,'(/,a)') 'Print transformation matrix (AOSO)'
   call print_r2mat(MatTr1,ntot,ntot,5,fileout,2)
   
   call get_aces2_matrix1(ntot,MatTr2,'alice_aces2_mat_soao.txt',fileout)
   write(fileout,'(/,a)') 'Print transformation matrix (SOAO)'
   call print_r2mat(MatTr2,ntot,ntot,5,fileout,2)


   if (.true.) then
      call xgemm('T','N',ntot,ntot,ntot,1.0D0,MatTr1,ntot,MatS,ntot,0.0D0,Mat1,ntot)
!     call xgemm('N','N',ntot,ntot,ntot,1.0D0,MatTr2,ntot,MatS,ntot,0.0D0,Mat1,ntot)
      call xgemm('N','N',ntot,ntot,ntot,1.0D0,Mat1,ntot,MatTr1,ntot,0.0D0,MatS,ntot)
      
      write(fileout,'(/,a)') 'Print Overlap (SOAO)'
!     call print_r2mat(Mat1,ntot,ntot,5,fileout,1)
      call print_r2mat(MatS,ntot,ntot,5,fileout,1)
      
!     do i = 1, ntot
!        do j = 1, ntot
!           MatS(i,j) = Mat1(i,j)
!        enddo
!     enddo
   endif

!  call xgemm('N','N',ntot,ntot,ntot,1.0D0,MatTr2,ntot,MatS,ntot,0.0D0,Mat1,ntot)
!  call xgemm('N','N',ntot,ntot,ntot,1.0D0,Mat1,ntot,MatTr2,ntot,0.0D0,MatS,ntot)
!  write(fileout,'(/,a)') 'Print Overlap (SOAO)'
!  call print_r2mat(MatS,ntot,ntot,5,fileout,1)
   

   if (.false.) then
   if (itype==1) then
      !normalize
      do j = 1, ntot
         temp = 0.0d0 
         do i = 1, ntot
            temp = temp + MatS(i,j)*MatS(i,j)
         enddo
         do i = 1, ntot
            MatS(i,j) = MatS(i,j)/temp
         enddo
      enddo

      write(fileout,'(/,a)') 'Print Overlap (normalized)'
      call print_r2mat(MatS,ntot,ntot,5,fileout,1)
   endif   
   endif   
end subroutine get_aces2_overlap


subroutine get_aces2_info(nirrep,occupya,numbasir,nlindep,molab,fileout)
   implicit none
   character(80) string
   integer nirrep,occupya(8),numbasir(8),nlindep(8)
   character(5) molab(8)
   integer i,fileinp,fileout

   fileinp = 21
   open(unit=fileinp,file="alice_aces2_info.txt",status='old',position='rewind')
   do while (.true.)
      read(fileinp,'(A80)',end=90) string
      if (string(1:7)=='*NIRREP') then 
         read(string(11:80),*) nirrep
         write(fileout,'(a,i5)') 'ACES2,NIRREP : ',nirrep
      elseif (string(1:8)=='*OCCUPYA') then 
         read(string(11:80),*) (occupya(i),i=1,nirrep)
         write(fileout,'(a,8i5)') 'ACES2,OCCUPYA : ',(occupya(i),i=1,nirrep)
      elseif (string(1:9)=='*NUMBASIR') then 
         read(string(11:80),*) (numbasir(i),i=1,nirrep)
         write(fileout,'(a,8i5)') 'ACES2,NUMBASIR : ',(numbasir(i),i=1,nirrep)
      elseif (string(1:6)=='*MOLAB') then 
         read(string(11:80),*) (molab(i),i=1,nirrep)
      elseif (string(1:8)=='*NLINDEP') then 
         read(string(11:80),*) (nlindep(i),i=1,nirrep)

      endif
   enddo
90 continue
   close(fileinp)
   return
end subroutine get_aces2_info

subroutine read_aces2_info(dim0,val,yesno,fileout,strinp)
   implicit none
   character(len=*) strinp
   character(80) string
   integer dim0,val(dim0)
   logical yesno
   integer i,fileinp,fileout

   fileinp=21
   open(unit=fileinp,file="alice_aces2_info.txt",status='old',position='rewind')

   yesno=.false.
   do while (.true.)
      read(fileinp,'(A80)',end=90) string

      if (string(2:7)==strinp .and. strinp=='NIRREP') then
         read(string(11:80),*) val(1)
         yesno=.true.
!        write(fileout,'(a,i5)') 'ACES2,NIRREP : ',val(1)

      elseif (string(2:8)==strinp .and. strinp=='OCCUPYA') then
         read(string(11:80),*) (val(i),i=1,dim0)
         yesno=.true.
!        write(fileout,'(a,8i5)') 'ACES2,OCCUPYA : ',(val(i),i=1,dim0)

      elseif (string(2:9)==strinp .and. strinp=='NUMBASIR') then
         read(string(11:80),*) (val(i),i=1,dim0)
         yesno=.true.
!        write(fileout,'(a,8i5)') 'ACES2,NUMBASIR : ',(val(i),i=1,dim0)

      elseif (string(2:8)==strinp .and. strinp=='NLINDEP') then
         read(string(11:80),*) (val(i),i=1,dim0)
         yesno=.true.

      elseif (string(2:9)==strinp .and. strinp=='NBASCART') then
         read(string(11:80),*) val(1)
         yesno=.true.
!        write(fileout,'(a,i5)') 'ACES2,NBASCART : ',val(1)

      elseif (string(2:8)==strinp .and. strinp=='NBASTOT') then
         read(string(11:80),*) val(1)
         yesno=.true.
!        write(fileout,'(a,i5)') 'ACES2,NBASTOT : ',val(1)

      endif

      if (yesno) goto 90
   enddo

90 continue
   close(fileinp)
   return
end subroutine read_aces2_info

subroutine get_aces2_basinfo(ntot,ao_nuc,ao_str,ao_ind,ao_ang, &
           fileread,fileout)
   implicit none
   character(len=*) fileread
   character(80) string
   integer ntot
   character(2), allocatable :: nuc(:)
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

!  do while (.true.)
!     read(fileinp,'(a80)',end=90) string
!     if (string(1:1)=='*') exit
!     read(string,*)
!  enddo

   i = 0
   do while (.true.) 
      read(fileinp,'(a80)',end=90) string
      if (string(1:1)=='*') exit      
      i = i + 1
      read(string,*) nuctype(i),nucnum(i,1),nuc(i),nucnum(i,2),nucnum(i,3), &
          geo(i,1),geo(i,2),geo(i,3) 
   enddo

   !AO list
   !write(fileout,'(a)') 'reading AO list'
   i = 0
   do while (.true.) 
      read(fileinp,'(a80)',end=90) string
      if (string(1:1)=='*') exit      
      !write(fileout,*) string
!                   |nucs nuc1 nuc2|ang1 ang2 angs|bas1 bas2 bas3 
      read(string,*) str1,ind1,ind2,ind3,ind4,str2,ind5,ind6,ind7 

!     do j = 1, ind6
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
!     enddo
   enddo
90 continue
   close(fileinp)

!  do i = 1, ntot-1
!     do j = i+1, ntot
!        if (ao_ind(j,2)/=1) then
!        enddo
!     enddo
!  enddo

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

end module mod_get_aces2
