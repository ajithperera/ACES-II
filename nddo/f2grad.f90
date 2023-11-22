subroutine f2grad(loop,fock,density,twoe)
use indices
use constants
implicit double precision(a-h,o-z)
double precision,dimension(:),intent(inout)::fock
double precision,dimension(:),intent(in)::density,twoe
double precision,dimension(:,:),allocatable::plocal,ploc
integer::mu,ka,la,nu
integer,dimension(2),intent(in)::loop
do mm=1,1
   i=loop(mm)
   istart=ifirst(i)
   iend=ilast(i)
   ni=iend-istart+1
   istart2=ifirst2(i)
   iend2=ilast2(i)

do nn=2,2
   j=loop(nn)
   jstart=ifirst(j)
   jend=ilast(j)
   nj=jend-jstart+1
   jstart2=ifirst2(j)
   jend2=ilast2(j)

  allocate(plocal(nj,ni))
      iplocal=0
   jplocal=0
   do itemp=jstart,jend
      iplocal=iplocal+1
      jplocal=0
   do jtemp=istart,iend
      jplocal=jplocal+1
      call pack1(jtemp,itemp,ij)
      plocal(iplocal,jplocal)=density(ij)
   end do
   end do
 
  jcount=0
   joffset=0
   ioffset=0
! loop over pairs on i and j and put in F(ka,mu) blocks 
do mu=istart,iend
   ioffset=ioffset+1
   joffset=0
!   if(ioffset>4)return
do ka=jstart,jend  ! so F(ka,mu) is the number we are after for this pair
total=zero
joffset=joffset+1
!if(joffset>4)return
!   the loop over mn is for the sum over nu lambda part so there are 16(max) of these for an sp basis set
      do la=1,ni
      do nu=1,nj
call pack2(jstart2+map(joffset,nu),istart2+map(ioffset,la),ij)
total=total+plocal(nu,la)*twoe(ij)
   end do
   end do
fock(mu+offset1(ka))=-half*total

end do
end do

deallocate(plocal)
end do
end do







end subroutine f2grad
