subroutine find
use indices
use constants
use tables
     use control
use scratch_array
implicit double precision (a-h,o-z)
interface

subroutine square(matrix,vector,idim1)
use indices
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(inout)::matrix
double precision,dimension(:),intent(in)::vector
integer,intent(in)::idim1
end subroutine square

 subroutine twoe_pair(i,j,t1)
 integer,intent(in)::i,j
 double precision,intent(inout)::t1
end subroutine twoe_pair

subroutine hcoregrad(loop)
integer,dimension(2),intent(in)::loop
end subroutine hcoregrad

subroutine f1grad(loop,h,fock,density,twoe)
double precision,dimension(:),intent(inout)::fock
integer,dimension(2),intent(in)::loop
double precision,dimension(:),intent(in)::density,twoe,h
end subroutine f1grad

subroutine f2grad(loop,fock,density,twoe)
integer,dimension(2),intent(in)::loop
double precision,dimension(:),intent(inout)::fock
double precision,dimension(:),intent(in)::density,twoe
end subroutine f2grad

subroutine trace(S,N,out)
double precision,dimension(:,:),intent(in)::S
integer::N
double precision,intent(out)::out
end subroutine trace

subroutine hcoregradsp(i)
integer,intent(in)::i
end subroutine hcoregradsp

      subroutine pack1(indi,indj,ij)
integer,intent(in)::indi,indj
integer,intent(inout)::ij
end subroutine pack1


end interface


double precision,allocatable,dimension(:,:)::grad
!double precision,intent(inout),dimension(:)::g
double precision::delta,xold,yold,zold
double precision,dimension(:,:),allocatable::pone,ptwo,hlocal  & ! i.e. p-local,h-local etc
,temp,temp2,temp3,pthr
double precision,dimension(:),allocatable::fock
double precision,dimension(2)::force
integer,dimension(2)::loop


!  matrix s3 holds the frozen density from the SCF calculation
allocate(grad(3,numat))
if(.not.save_tree)then
print*,'Computing gradient using RHF determinant...'
end if
call cpusec(time1)
delta=1d-3 ! smaller values of delta for the finite diff cause num. inaccuracies
eold=etotal

allocate(fock(ndim1))
allocate(scrvec(ndim1))
allocate(scrmat(num_basis,num_basis))
grad=zero
fock=zero
do i=1,numat-1

i1=ifirst(i)
i2=ilast(i)
aa=zero
bb=zero
cc=zero
dd=zero
 

 do j=i+1,numat
j1=ifirst(j)
j2=ilast(j)
ni=nbas(species(i))
nj=nbas(species(j))
loop(1)=i
loop(2)=j

xold=x(i)
x(i)=xold+delta

!vdbh=zero
call twoe_pair(i,j,t1)
call hcoregrad(loop)
!vdbfock=zero
call f1grad(loop,h,fock,s3,twoe)
call  f2grad(loop,fock,s3,twoe)


energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*(two*h(ij)+fock(ij))
   end do
   end do

   do itemp=j1,j2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*(two*h(ij)+fock(ij))
   end do
   end do

   do itemp=i1,i2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
      energy=energy+two*s3(ij)*(two*h(ij)+fock(ij))
   end do
   end do



force(1)=half*energy+t1



x(i)=xold-delta
!vdbh=zero
call twoe_pair(i,j,t1)
call hcoregrad(loop)
!vdbfock=zero
call f1grad(loop,h,fock,s3,twoe)
call  f2grad(loop,fock,s3,twoe)


energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*(two*h(ij)+fock(ij))
   end do
   end do

   do itemp=j1,j2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*(two*h(ij)+fock(ij))
   end do
   end do

   do itemp=i1,i2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
      energy=energy+two*s3(ij)*(two*h(ij)+fock(ij))
   end do
   end do

force(2)=half*energy+t1


der=(force(1)-force(2))/(two*delta)
grad(1,i)=grad(1,i)+der
grad(1,j)=grad(1,j)-der

x(i)=xold


yold=y(i)
y(i)=yold+delta

!vdbh=zero
call twoe_pair(i,j,t1)
call hcoregrad(loop)
!vdbfock=zero
call f1grad(loop,h,fock,s3,twoe)
call  f2grad(loop,fock,s3,twoe)


energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*(two*h(ij)+fock(ij))
  end do
   end do

   do itemp=j1,j2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*(two*h(ij)+fock(ij))
   end do
   end do

   do itemp=i1,i2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
      energy=energy+two*s3(ij)*(two*h(ij)+fock(ij))
   end do
   end do



force(1)=half*energy+t1



y(i)=yold-delta
!vdbh=zero
call twoe_pair(i,j,t1)
call hcoregrad(loop)
!vdbfock=zero
call f1grad(loop,h,fock,s3,twoe)
call  f2grad(loop,fock,s3,twoe)


energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*(two*h(ij)+fock(ij))
   end do
   end do

   do itemp=j1,j2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*(two*h(ij)+fock(ij))
   end do
   end do

   do itemp=i1,i2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
      energy=energy+two*s3(ij)*(two*h(ij)+fock(ij))
   end do
   end do

force(2)=half*energy+t1


der=(force(1)-force(2))/(two*delta)
grad(2,i)=grad(2,i)+der
grad(2,j)=grad(2,j)-der

y(i)=yold


zold=z(i)
z(i)=zold+delta

!vdbh=zero
call twoe_pair(i,j,t1)
call hcoregrad(loop)
!vdbfock=zero
call f1grad(loop,h,fock,s3,twoe)
call  f2grad(loop,fock,s3,twoe)


energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*(two*h(ij)+fock(ij))
   end do
   end do

   do itemp=j1,j2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*(two*h(ij)+fock(ij))
   end do
   end do

   do itemp=i1,i2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
      energy=energy+two*s3(ij)*(two*h(ij)+fock(ij))
   end do
   end do



force(1)=half*energy+t1



z(i)=zold-delta
!vdbh=zero
call twoe_pair(i,j,t1)
call hcoregrad(loop)
!vdbfock=zero
call f1grad(loop,h,fock,s3,twoe)
call  f2grad(loop,fock,s3,twoe)


energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*(two*h(ij)+fock(ij))
   end do
   end do

   do itemp=j1,j2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*(two*h(ij)+fock(ij))
   end do
   end do

   do itemp=i1,i2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
      energy=energy+two*s3(ij)*(two*h(ij)+fock(ij))
   end do
   end do

force(2)=half*energy+t1


der=(force(1)-force(2))/(two*delta)
grad(3,i)=grad(3,i)+der
grad(3,j)=grad(3,j)-der

z(i)=zold







end do




end do



if(sparkles)then
do i=1,numat
i1=ifirst(i)
i2=ilast(i)

xold=x(i)
x(i)=xold+delta
enuc=zero
h=zero
call hcoregradsp(i)

energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*two*h(ij)
   end do
   end do




force(1)=half*energy+enuc








x(i)=xold-delta
enuc=zero
h=zero
call hcoregradsp(i)
energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*two*h(ij)
   end do
   end do





force(2)=half*energy+enuc

der=(force(1)-force(2))/(two*delta)
grad(1,i)=grad(1,i)+der


x(i)=xold

yold=y(i)
y(i)=yold+delta
enuc=zero
h=zero
call hcoregradsp(i)
energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*two*h(ij)
   end do
   end do





force(1)=half*energy+enuc




y(i)=yold-delta
enuc=zero
h=zero
call hcoregradsp(i)
energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*two*h(ij)
   end do
   end do




force(2)=half*energy+enuc



der=(force(1)-force(2))/(two*delta)
grad(2,i)=grad(2,i)+der


y(i)=yold


zold=z(i)
z(i)=zold+delta
enuc=zero
h=zero
call hcoregradsp(i)
energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*two*h(ij)
   end do
   end do




force(1)=half*energy+enuc
z(i)=zold-delta
enuc=zero
h=zero
call hcoregradsp(i)
energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*two*h(ij)
   end do
   end do



force(2)=half*energy+enuc



der=(force(1)-force(2))/(two*delta)
grad(3,i)=grad(3,i)+der


z(i)=zold

end do


end if











grad=grad*627.51d0/27.21d0


if(allocated(g))then
do i=1,numat
      ISTART=3*(I-1)+1
      g(istart)=grad(1,i)
      g(istart+1)=grad(2,i)
      g(istart+2)=grad(3,i)
end do
end if
















rnorm=zero
if(.not.save_tree)then
print*,'                     Gradient'
print*,'--------------------------------------------------------------'
do i=1,numat
write(*,10)i,grad(1,i),'     ',grad(2,i),'     ',grad(3,i)
rnorm=rnorm+grad(1,i)**2+grad(2,i)**2+grad(3,i)**2
end do 

rnorm=dsqrt(rnorm)
print*,'||F|| =',rnorm
print*,'--------------------------------------------------------------'
10 format(i4,f20.12,f20.12,f20.12)
call cpusec(time2)
write(*,202)time2-time1
202 format(' Elapsed time for gradient calculation (sec) =  ',F10.3)
end if

deallocate(grad)
deallocate(scrvec)
deallocate(scrmat)

end subroutine find



