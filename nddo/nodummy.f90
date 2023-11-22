subroutine nodummy
use tables
use constants
double precision,dimension(:),allocatable::oldx,oldy,oldz,oldxsp,oldysp,oldzsp
integer,dimension(:),allocatable::oldzeff,oldzeffsp

allocate(oldx(ninput))
allocate(oldy(ninput))
allocate(oldz(ninput))
allocate(oldzeff(ninput))
numat=0
nsparkle=0
irow=0

if(allocated(xsparkle))deallocate(xsparkle)
if(allocated(ysparkle))deallocate(ysparkle)
if(allocated(zsparkle))deallocate(zsparkle)
if(allocated(zeffsp))deallocate(zeffsp)
allocate(xsparkle(1))
allocate(ysparkle(1))
allocate(zsparkle(1))
allocate(zeffsp(1))

do i=1,ninput

 if(zstore(i)>0  .and. zstore(i).le.85)then

 numat=numat+1
 irow=irow+1
 oldx(irow)=x(i)
 oldy(irow)=y(i)
 oldz(irow)=z(i)
 oldzeff(irow)=zstore(i)

elseif(zstore(i).gt.85)then ! pluck out the sparkle atoms
nsparkle=nsparkle+1
allocate(oldxsp(nsparkle))
allocate(oldysp(nsparkle))
allocate(oldzsp(nsparkle))
allocate(oldzeffsp(nsparkle))
oldxsp=xsparkle
oldysp=ysparkle
oldzsp=zsparkle
oldzeffsp=zeffsp
deallocate(xsparkle)
deallocate(ysparkle)
deallocate(zsparkle)
deallocate(zeffsp)
allocate(xsparkle(nsparkle))
allocate(ysparkle(nsparkle))
allocate(zsparkle(nsparkle))
allocate(zeffsp(nsparkle))
xsparkle=oldxsp
ysparkle=oldysp
zsparkle=oldzsp
zeffsp=oldzeffsp
xsparkle(nsparkle)=x(i)
ysparkle(nsparkle)=y(i)
zsparkle(nsparkle)=z(i)
zeffsp(nsparkle)=zstore(i)
deallocate(oldxsp)
deallocate(oldysp)
deallocate(oldzsp)
deallocate(oldzeffsp)
 end if
end do

deallocate(x)
deallocate(y)
deallocate(z)
deallocate(zeff)

allocate(x(numat))
allocate(y(numat))
allocate(z(numat))
allocate(zeff(numat))
x=oldx
y=oldy
z=oldz
zeff=oldzeff
deallocate(oldx)
deallocate(oldy)
deallocate(oldz)
deallocate(oldzeff)


zeffsp=zeffsp-85


end subroutine nodummy
