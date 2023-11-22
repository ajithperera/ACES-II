subroutine f1uhf(fock,density,abdens)
use indices
use constants
use tables
implicit double precision(a-h,o-z)
double precision,dimension(:),intent(inout)::fock
double precision,dimension(:),intent(in)::density,abdens
double precision,dimension(:,:),allocatable::plocal,ploc,pabloc
fock=h
tot=zero
do i=1,numat
   istart=ifirst(i)
   iend=ilast(i)
   ni=iend-istart+1
   istart2=ifirst2(i)
   iend2=ilast2(i)
   allocate(plocal(ni,ni))
   allocate(pabloc(ni,ni))
      iplocal=0
   jplocal=0
   do itemp=istart,iend
      iplocal=iplocal+1
      jplocal=0
   do jtemp=istart,iend
      jplocal=jplocal+1
      call pack1(jtemp,itemp,ij)
      plocal(iplocal,jplocal)=density(ij)
      pabloc(iplocal,jplocal)=abdens(ij)
   end do
   end do




  ! j is index for loop over basis fn's on atom i
   icount=0
   mm=0
    ispot=0
    jspot=0
    
   do j=istart,iend
      
      total=zero
      icount=icount+1
          do k=1,ni
          do l=1,ni

call pack2(istart2+map(icount,icount),istart2+map(l,k),ij)
call pack2(istart2+map(icount,l),istart2+map(icount,k),kl)
 total=total+plocal(l,k)*twoe(ij)-pabloc(l,k)*twoe(kl) 

end do
end do

fock(j+offset1(j))=fock(j+offset1(j))+total
   ! now, loop over uu terms on center B

             total=zero
           do k=1,numat
           if(k==i)goto 10
           kstart=ifirst(k)
           klast=ilast(k)
            nk=klast-kstart+1
           kstart2=ifirst2(k)
           allocate(ploc(nk,nk))
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
              
                  call pack2(kstart2+map(m,l),icol,ij)
                  total=total+ploc(m,l)*twoe(ij)



               end do
               end do

              deallocate(ploc) 
10         continue
          
           end do
fock(j+offset1(j))=fock(j+offset1(j))+total
end do
deallocate(plocal)
deallocate(pabloc)
! this ends the F(u,u) section
! do the f(u,v) section
! loop over the uv pairs
 ! hydrogen atoms screw this up so check for that special case
  if(ni>1)then
allocate(ploc(ni,ni))
allocate(pabloc(ni,ni))
   iplocal=0
   jplocal=0
   do itemp=istart,iend
      iplocal=iplocal+1
      jplocal=0
   do jtemp=istart,iend
      jplocal=jplocal+1
      call pack1(jtemp,itemp,ij)
      ploc(iplocal,jplocal)=density(ij)
      pabloc(iplocal,jplocal)=abdens(ij)
   end do
   end do

! these loop indices pick up the offdiagonal one center blocks
  ispot=0
    do l=istart,iend-1
    ispot=ispot+1
    mspot=ispot+1
   do m=l+1,iend
total=zero
do k=1,ni
do n=1,ni
call pack2(istart2+map(ispot,mspot),istart2+map(n,k),ij)
call pack2(istart2+map(ispot,n),istart2+map(mspot,k),kl)
    total=total+ploc(n,k)*twoe(ij)-pabloc(n,k)*twoe(kl) 


end do
end do




      do k=1,numat
      if(k==i)goto 90

           kstart=ifirst(k)
           klast=ilast(k)
           nk=klast-kstart+1
           kstart2=ifirst2(k)
           allocate(plocal(nk,nk))

           kcount=1
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




           do ll=1,nk
           do mm=1,nk
              call pack2(kstart2+map(mm,ll),istart2+map(mspot,ispot),kl)
total=total+plocal(mm,ll)*twoe(kl)

               kcount=kcount+1
               end do
               end do




               deallocate(plocal)
      90 continue
      end do
       
      fock(l+offset1(m))=fock(l+offset1(m))+total
mspot=mspot+1
    end do
   end do
deallocate(ploc)
deallocate(pabloc)
end if

     


end do
 
! this end the loop over atoms

end subroutine f1uhf
