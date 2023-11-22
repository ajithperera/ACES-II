subroutine f1grad(loop,h,fock,density,twoe)
use indices
use constants
implicit double precision(a-h,o-z)
double precision,dimension(:),intent(inout)::fock
double precision,dimension(:),intent(in)::density,twoe,h
!double precision,dimension(4,4)::plocal,ploc
double precision,dimension(:,:),allocatable::plocal,ploc
integer,dimension(2),intent(in)::loop
!fock=h
tot=zero
do m1=1,2
i=loop(m1)
!   ni=nbas(i)
   istart=ifirst(i)
   iend=ilast(i)
   ni=iend-istart+1
   istart2=ifirst2(i)
   iend2=ilast2(i)
   icount=0
   mm=0
    ispot=0
    jspot=0
   do j=istart,iend
      icount=icount+1
      total=zero
   ! now, loop over uv terms on center B

             total=zero
           do nn=1,2
            k=loop(nn)  
           if(k==i)goto 10
!           nk=nbas(k)
           kstart=ifirst(k)
           klast=ilast(k)
            nk=klast-kstart+1
           kstart2=ifirst2(k)
           allocate(ploc(nk,nk))
!p           ploc(1:nk,1:nk)=density(kstart:klast,kstart:klast)
   iplocal=0
   jplocal=0
   do itemp=kstart,klast
      iplocal=iplocal+1
      jplocal=0
   do jtemp=kstart,klast
      jplocal=jplocal+1
      call pack1(jtemp,itemp,ij)
      ploc(iplocal,jplocal)=density(ij)
   end do
   end do

           kcount=1
           icol=istart2+map(icount,icount)
               do l=1,nk
               do m=1,nk
           
              !total=total+ploc(m,l)*twoe(kstart2+map(m,l),icol)
               call pack2(kstart2+map(m,l),icol,ij)
                  total=total+ploc(m,l)*twoe(ij) 
               end do
               end do

               deallocate(ploc) 
10         continue
          
           end do

!fock(j,j)=total
fock(j+offset1(j))=total
           
end do

! this ends the F(u,u) section
! do the f(u,v) section
! loop over the uv pairs
 ! hydrogen atoms screw this up so check for that special case
  if(ni>1)then
  ispot=0
    do l=istart,iend-1
    ispot=ispot+1
    mspot=ispot+1
   do m=l+1,iend
    total=zero
      do l2=1,2
         k=loop(l2)
      if(k==i)goto 90

           kstart=ifirst(k)
           klast=ilast(k)
           nk=klast-kstart+1
           kstart2=ifirst2(k)
           allocate(plocal(nk,nk))
!p           plocal(1:nk,1:nk)=density(kstart:klast,kstart:klast)
           iplocal=0
   jplocal=0
   do itemp=kstart,klast
      iplocal=iplocal+1
      jplocal=0
   do jtemp=kstart,klast
      jplocal=jplocal+1
      call pack1(jtemp,itemp,ij)
      plocal(iplocal,jplocal)=density(ij)
   end do
   end do

           kcount=1
           do ll=1,nk
           do mm=1,nk
 call pack2(kstart2+map(mm,ll),istart2+map(mspot,ispot),kl)
total=total+plocal(mm,ll)*twoe(kl)               
               end do
               end do




deallocate(plocal)
      90 continue
      end do
fock(l+offset1(m))=total
mspot=mspot+1
    end do
   end do
end if

     


end do
 
! this end the loop over atoms
end subroutine f1grad
