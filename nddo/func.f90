subroutine func(x,value,L,condon,D)
use constants
implicit double precision (a-h,o-z)
double precision,intent(in)::x,condon,D
double precision,intent(inout)::value
integer,intent(in)::L

if(L==1)then
xx=x+x
value=half * ( one / ( xx )  - one / sqrt ( four * D*D +  xx * xx )  ) - condon
elseif(L==2)then
xx=x+x
xxxx=xx*xx
value = pt25 /  sqrt(xxxx)   +  pt25 / sqrt ( 8.0d0 * D * D + xxxx ) - half / sqrt(four * D * D + xxxx) - condon

end if
end subroutine func
