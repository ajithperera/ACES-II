subroutine scfopt(x0,y0,z0,gout,natoms,ff)
use tables
use control
use indices
use constants
implicit double precision(a-h,o-z)
interface

subroutine find
end subroutine find

subroutine finduhf
end subroutine finduhf
end interface

double precision,intent(inout) ::ff
double precision, dimension(natoms),intent(inout) ::x0,y0,z0
 double precision, dimension(3*natoms),intent(inout) ::gout
!size(x0)
x=x0
y=y0
z=z0
if(optimize)then
write(*,*)' ATOMIC NO.','      X   ','        Y   ','         Z'
do i=1,numat
write(*,300)zeff(i),x(i),y(i),z(i)
end do
300 format(i8,f13.5,f13.5,f13.5)
end if
deallocate(ifirst)
deallocate(ifirst2)
deallocate(ilast)
deallocate(ilast2)
call two_electron
call hcore

TRIAL=.false.
keep=.true.
if(restricted)then
call rhf
else
call uhf
end if
TRIAL=.false.
keep=.false.

if(restricted)then
call find
else
call finduhf
end if

gout=g
ff=etotal

return 
end


