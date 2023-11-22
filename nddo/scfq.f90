subroutine scfq(q0,natoms,ff,stat1,stat2)
use tables
use control
use indices
use constants
implicit double precision(a-h,o-z)
interface

subroutine find
end subroutine find

subroutine geometry
end subroutine geometry
end interface
integer,intent(in)::natoms
double precision,intent(inout) ::ff
! double precision, dimension(3*natoms),intent(inout) ::gout
double precision,dimension(natoms,3),intent(in)::q0
logical::stat1,stat2
!size(x0)
q=q0
if(stat2)then
write(*,*)'GEOMETRY AT CURRENT ITERATE'
     do i=1,ninput
     write(*,301)zstore(i),q(i,1),opt(i,1),q(i,2)*180.0d0/pi,opt(i,2),q(i,3)*180.0d0/pi,opt(i,3),ref(i,1),ref(i,2),ref(i,3)
     end do
301 format(i3,f10.5,i3,f10.5,i3,f10.5,i3,i3,i3,i3)
end if
!TRIAL=status
!keep=stat2
call geometry
call nodummy




deallocate(ifirst)
deallocate(ifirst2)
deallocate(ilast)
deallocate(ilast2)
call two_electron
call hcore

TRIAL=stat1
keep=stat2
if(restricted)then
call rhf
else
call uhf
end if
!TRIAL=.false.
!keep=.false.

!call find
!gout=g
ff=etotal

!call find
return 
end


