module mod_aovec
contains

subroutine aovec_trans_sym(nsph,trmat,fileout)
   use mod_print
   implicit none
   integer nsph,fileout,fileinp,i,num(5),ncart
   double precision trmat(nsph,nsph),val
   character(80) string

   fileinp=52
   trmat(1:nsph,1:nsph) = 0.0d0
   open(unit=fileinp,file='ACES2_BASINFO2',status='old',position='rewind')

   write(fileout,'(a)') 'read ACES2_BASINFO2'
   do while(.true.)
      read(fileinp,'(a80)',end=10) string
      !num: 1.No(cart)/2(atom)/3(l)/4(lz)/5.No(sph)
      read(string(1:35),*) (num(i),i=1,5),val,ncart

      if (num(3)<=nsph) then 
         !trmat(sph(sym),sph(org))
         trmat(num(5),num(3)) = val
         write(fileout,'(2i5,a,f5.1)') num(5),num(3),' = ',val
         !read(string(36:80),*) (sphorb(i),i=1,ncart)
      endif
   enddo
10 continue
   close(fileinp)

   write(fileout,'(/,a)') '* trmat (sph->sph)'
   call print_rmat(trmat,nsph*nsph,nsph,nsph,fileout,2)
end subroutine aovec_trans_sym


subroutine aovec_sph2cart_addfn(ncart,nsph,aoind,aoind2,fileout)
   use mod_fun
   implicit none
   integer ncart,nsph,aoind(nsph,5),aoind2(ncart,6),fileout
   integer i,j,icount
   icount=0
   do i = 1, nsph
      icount=icount+1
      do j=1,5
         aoind2(icount,j)=aoind(i,j)
         aoind2(icount,6)=icount
      enddo

      if ((aoind(i,3)==3) .and. (aoind(i,4)==4)) then
         !d ordering 2/5/1/3/4 -> d(-2)/d(-1)/d(0)/d(1)/d(2)
         !add d(6) orbital after d(4)
         icount=icount+1
         aoind2(icount,1)=0
         aoind2(icount,2)=aoind(i,2)
         aoind2(icount,3)=aoind(i,3)
         aoind2(icount,4)=6
         aoind2(icount,5)=aoind(i,5)
         aoind2(icount,6)=icount
      endif
   enddo

   write(fileout,'(/,a)') '* aolab (add functions)'
   do i = 1,ncart
      write(fileout,'(6i5)') (aoind2(i,j),j=1,6)
   enddo
end subroutine aovec_sph2cart_addfn


subroutine aovec_make_trmat(nsph,ncart,aoind2,trmat,fileout)
   use mod_fun
   use mod_print
   implicit none
   integer ncart,nsph,aoind2(ncart,6),fileout
   double precision trmat(nsph,ncart)
   integer i,j,i2,ind

   trmat(1:nsph,1:ncart)=0.0d0

   do i = 1, ncart
      if (aoind2(i,3)<=2) then
         i2=aoind2(i,6)
         trmat(aoind2(i,1),i2) = 1.0d0
      elseif (aoind2(i,3)==3) then
         !d ordering 2/5/1/3/4/6 -> d(-2)/d(-1)/d(0)/d(1)/d(2)/d(xx)
!           1      2    3    4      5    6
!           2      5    1    3      4    6
!           d0     d-2  d1   d2     d-1  xx
!        / '-1/6', '0', '0',' 1/2', '0', '0', &  !1
!          '   0', '1', '0','   0', '0', '0', &  !2
!          '   0', '0', '1','   0', '0', '0', &  !3
!          '-1/6', '0', '0','-1/2', '0', '0', &  !4
!          '   0', '0', '0','   0', '1', '0', &  !5
!          ' 1/3', '0', '0','   0', '0', '0' /   !6

         i2=aoind2(i,6)
         ind=aoind2(i,4)
         if (ind==1) then !dxx
!           trmat(aoind2(i-2,1),i2) = -1.0/sqrt(6.0)
!           trmat(aoind2(i+1,1),i2) =  1.0/sqrt(2.0)
!           trmat(aoind2(i-2,1),i2) = -1.0/6.0
!           trmat(aoind2(i+1,1),i2) =  1.0/2.0
            trmat(aoind2(i  ,1),i2) = -1.0/4.0
            trmat(aoind2(i+3,1),i2) =  3.0/4.0
!           trmat(aoind2(i  ,1),i2) = -1.0/6.0
!           trmat(aoind2(i+3,1),i2) =  1.0/2.0
         elseif (ind==2) then !dxy
!           trmat(aoind2(i+1,1),i2) =  1.0
            trmat(aoind2(i  ,1),i2) =  1.0
         elseif (ind==3) then !dxz
!           trmat(aoind2(i-1,1),i2) =  1.0
            trmat(aoind2(i  ,1),i2) =  1.0
         elseif (ind==4) then !dyy
!           trmat(aoind2(i-4,1),i2) = -1.0/sqrt(6.0)
!           trmat(aoind2(i-1,1),i2) = -1.0/sqrt(2.0)
!           trmat(aoind2(i-4,1),i2) = -1.0/6.0
!           trmat(aoind2(i-1,1),i2) = -1.0/2.0
            trmat(aoind2(i-3,1),i2) = -1.0/4.0
            trmat(aoind2(i  ,1),i2) = -3.0/4.0
!           trmat(aoind2(i-3,1),i2) = -1.0/6.0
!           trmat(aoind2(i  ,1),i2) = -1.0/2.0
         elseif (ind==5) then !dyz
!           trmat(aoind2(i+3,1),i2) =  1.0
            trmat(aoind2(i  ,1),i2) =  1.0
         elseif (ind==6) then !dzz
            trmat(aoind2(i-5,1),i2) =  1.0
!           trmat(aoind2(i-5,1),i2) =  1.0/3.0
         endif
      endif
   enddo

   write(fileout,'(/,a)') '* trmat'
   call print_rmat(trmat,ncart*nsph,nsph,ncart,fileout,2)
end subroutine aovec_make_trmat


subroutine aovec_aosph_nw2ac(ntot,natom,aoind,aostr,geoind,trmat,fileout)
   use mod_fun
   use mod_print
   implicit none
   integer ntot,natom,fileout
   integer aoind(ntot,5),geoind(natom,3)
   double precision trmat(ntot,ntot)
   character(4) aostr(ntot,2)
   integer aotmp(5)
   integer i,j,k,b1,b2,c1,c2,k1,k2,l1,l2,m1,m2

!  aoind(i,j,k,l,m) i:nw-order/j:atom/k:ang/l:ang(lz)/m:ang(lz,ind)
!  geoind(a,b,c)    a:atom/b:atom-group/c:group-ind
!  aces-ind         b->k->l->m->c
   do i = 1,ntot-1
   do j = i+1,ntot
      ! atom group
      b1=geoind(aoind(i,2),2)
      b2=geoind(aoind(j,2),2)
      if (b2.lt.b1) then
         call rot_aoind(ntot,aoind,i,j)
      elseif (b2.eq.b1) then
         ! angular momentum
         k1=aoind(i,3)
         k2=aoind(j,3)
         if (k2.lt.k1) then
            call rot_aoind(ntot,aoind,i,j)
         elseif (k2.eq.k1) then
            ! angular momentum (z-axis)
            l1=aoind(i,4)
            l2=aoind(j,4)
            if (l2.lt.l1) then
               call rot_aoind(ntot,aoind,i,j)
            elseif (l2.eq.l1) then
               ! angular momentum (z-axis,index)
               m1=aoind(i,5)
               m2=aoind(j,5)
               if (m2.lt.m1) then
                  call rot_aoind(ntot,aoind,i,j)
               elseif (m2.eq.m1) then
                  ! atom group index
                  c1=aoind(i,2)
                  c2=aoind(j,2)
                  if (c2.lt.c1) then
                     call rot_aoind(ntot,aoind,i,j)
                  endif
               endif
            endif
         endif
      endif
   enddo
   enddo

   write(fileout,'(/,a)') '* aolab (nwchem rotated)'
   do i = 1,ntot
      write(fileout,'(5i5,2a4)') (aoind(i,j),j=1,5), &
           (aostr(aoind(i,1),k),k=1,2)
   enddo

   trmat(1:ntot,1:ntot) = 0.0d0
   do i = 1,ntot
      trmat(aoind(i,1),i) = 1.0d0  
   enddo

   write(fileout,'(/,a)') '* trmat (nwchem rotated)'
   call print_rmat(trmat,ntot*ntot,ntot,ntot,fileout,2)
end subroutine aovec_aosph_nw2ac

subroutine aovec_trans_total(dim1,dim2,dim3,trmat1,trmat2,trmat)
   use mod_print
   implicit none
   integer dim1,dim2,dim3,i,j,k
   double precision trmat1(dim1,dim3)
   double precision trmat2(dim2,dim3)
   double precision trmat(dim1,dim2)

   trmat(1:dim1,1:dim2) = 0.0d0
   do i = 1, dim1 
   do j = 1, dim2 
   do k = 1, dim3 
      trmat(i,j) = trmat(i,j) + trmat1(i,k)*trmat2(j,k) 
   enddo
   enddo
   enddo
end subroutine aovec_trans_total

subroutine aovec_rotate_sph(nsph,ncart,movec,trmat,fileout)
!  (1) trmat(nsph,ncart)' x movec(nsph,nsph)   -> movec2(ncart,nsph)
!  (2) AO2SO(nsph,ncart)  x movec2(ncart,nsph) -> movec(nsph,nsph)
   use mod_get_aces2
   use mod_print
   implicit none
   integer nsph,ncart,fileout
   double precision trmat(nsph,ncart),movec(nsph,nsph)
   integer i,j,k
   double precision, allocatable :: movec2(:,:)

   allocate(movec2(ncart,nsph)) 

!  (1) trmat(nsph,ncart)' x movec(nsph,nsph)   -> movec2(ncart,nsph)
   movec2(1:ncart,1:nsph)=0.0d0
   do i = 1, ncart
      do j = 1, nsph
         do k = 1, nsph
            movec2(i,j) = movec2(i,j) + trmat(k,i)*movec(k,j)
         enddo
      enddo
   enddo

!  read AO2SO ('alice_aces2_mat_ao2so.txt')
   call get_aces2_matrix2(nsph,ncart,trmat,'alice_aces2_mat_ao2so.txt',fileout)

!  (2) AO2SO(nsph,ncart)  x movec2(ncart,nsph) -> movec(nsph,nsph)
   movec(1:nsph,1:nsph)=0.0d0
   do i = 1, nsph
      do j = 1, nsph
         do k = 1, ncart
            movec(i,j) = movec(i,j) + trmat(i,k)*movec2(k,j)
         enddo
      enddo
   enddo

   deallocate(movec2) 

   write(fileout,'(/,a)') '* movec (nwchem->aces2)'
   call print_rmat(movec,nsph*nsph,nsph,nsph,fileout,2)
end subroutine aovec_rotate_sph

end module mod_aovec
