subroutine trace(S,N,out)
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(in)::S
integer::N
double precision,intent(out)::out
out=0.0d0
do i=1,N
out=out+S(i,i)
end do
return
end
