subroutine factorial(i,total)
implicit double precision(a-h,o-z)
integer,intent(inout)::total
integer,intent(in)::i
total=1
do k=i,2,-1
total=total*k
end do
return
end subroutine factorial
