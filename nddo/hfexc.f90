subroutine hfexc(coulomb,exchange,ekinetic)
use scratch_array
use constants
use tables
use indices

implicit double precision (a-h,o-z)


interface
subroutine trace(S,N,out)
double precision,dimension(:,:),intent(in)::S
integer::N
double precision,intent(out)::out
end subroutine trace

subroutine square(matrix,vector,idim1)
use indices
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(inout)::matrix
double precision,dimension(:),intent(in)::vector
integer,intent(in)::idim1
end subroutine square
end interface

double precision,intent(inout)::coulomb,exchange,ekinetic
double precision,dimension(:,:),allocatable::pa,pb,pab
integer::si
coulomb=zero
hfx=zero
do i=1,numat
ni=nbas(species(i))
istart=ifirst(i)
iend=ilast(i)
istart2=ifirst2(i)
iend2=ilast2(i)
allocate(pa(ni,ni))
!pa(1:ni,1:ni)=s3(istart:iend,istart:iend)
   iplocal=0
   jplocal=0
   do itemp=istart,iend
      iplocal=iplocal+1
      jplocal=0
   do jtemp=istart,iend
      jplocal=jplocal+1
      call pack1(jtemp,itemp,ij)
      pa(iplocal,jplocal)=s3(ij)
   end do
   end do
do j=1,numat
nj=nbas(species(j))
jstart=ifirst(j)
jend=ilast(j)
jstart2=ifirst2(j)
jend2=ilast2(j)
allocate(pb(nj,nj))

allocate(pab(ni,nj))


!pb(1:nj,1:nj)=s3(jstart:jend,jstart:jend)
   iplocal=0
   jplocal=0
   do itemp=jstart,jend
      iplocal=iplocal+1
      jplocal=0
   do jtemp=jstart,jend
      jplocal=jplocal+1
      call pack1(jtemp,itemp,ij)
      pb(iplocal,jplocal)=s3(ij)
   end do
   end do
!pab(1:ni,1:nj)=s3(istart:iend,jstart:jend)
      iplocal=0
   jplocal=0
   do itemp=istart,iend
      iplocal=iplocal+1
      jplocal=0
   do jtemp=jstart,jend
      jplocal=jplocal+1
      call pack1(jtemp,itemp,ij)
      pab(iplocal,jplocal)=s3(ij)
   end do
   end do

do mu=1,ni
do nu=1,ni
 do la=1,nj
 do si=1,nj
!value=twoe(jstart2+map(si,la),istart2+map(nu,mu))
call pack2(jstart2+map(si,la),istart2+map(nu,mu),ij)
value=twoe(ij)
coulomb=coulomb+pa(nu,mu)*pb(si,la)*value
hfx=hfx+pab(mu,la)*pab(nu,si)*value
end do
end do
end do
end do

deallocate(pb)
deallocate(pab)

end do
deallocate(pa)
end do
coulomb=coulomb/two
hfx=hfx/four
allocate(scrmat(num_basis,num_basis))
allocate(scrvec(ndim1))
scrvec=s3*(H)
call square(scrmat,scrvec,num_basis)
ekinetic=sum(scrmat)
!call xgemm( 'N', 'N', num_basis, num_basis,num_basis, one,s3 &
! , num_basis,H , num_basis,zero,scrmat,num_basis )
!call trace(scrmat,num_basis,ekinetic)
exchange=-hfx
deallocate(scrmat)
deallocate(scrvec)
end subroutine hfexc
