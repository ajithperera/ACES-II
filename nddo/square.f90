subroutine square(matrix,vector,idim1)
use indices
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(inout)::matrix
double precision,dimension(:),intent(in)::vector
integer,intent(in)::idim1
! idim1 is dimension of square matrix
! idim2 is length of vector
jindex=0
do i=1,idim1
do j=1,i
jindex=jindex+1
matrix(j,i)=vector(jindex)
matrix(i,j)=matrix(j,i)
end do
end do


!do i=1,idim1-1
!do j=i+1,idim1
!matrix(j,i)=matrix(i,j)
!end do
!end do
end subroutine square
