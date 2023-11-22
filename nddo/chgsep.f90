subroutine chgsep(L,out,condon,D)
use constants
implicit double precision(a-h,o-z)
interface

subroutine limit(low,high,L,condon,D)
double precision,intent(inout)::low,high
double precision,intent(in)::condon,D
integer,intent(in)::L
end subroutine limit
end interface
integer,intent(in)::L
double precision,intent(inout)::out
double precision,intent(in)::condon,D
double precision::low

if(L==0)then
out=half/condon
else
call limit(low,high,L,condon,D)
call rtbis(out,low,high,1D-10,L,condon,D)
end if


end subroutine chgsep




