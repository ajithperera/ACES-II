subroutine limit(low,high,L,condon,D)
use constants
implicit double precision(a-h,o-z)
interface
subroutine func(x,value,L,condon,D)
double precision,intent(in)::x,condon,D
double precision,intent(inout)::value
integer,intent(in)::L
end subroutine func
end interface


double precision,intent(inout)::low,high
double precision,intent(in)::condon,D
integer,intent(in)::L
delta=1d-4
a=0.0001d0
b=a + delta
call func(a,aout,L,condon,D)
call func(b,bout,L,condon,D)
if(aout>zero .and. bout>zero)goto 20

do
delta=delta/10.0d0
a=delta
b=a+delta
call func(a,aout,L,condon,D)
call func(b,bout,L,condon,D)
if(aout>zero .and. bout>zero)exit
end do






20 continue
do 
  call func(a,aout,L,condon,D)
  call func(b,bout,L,condon,D)
  if ( (aout>zero) .and. (bout<zero) )then
     
  low=a
  high=b
  return
  end if
  a=a+delta
  b=a+delta
end do

end subroutine limit



