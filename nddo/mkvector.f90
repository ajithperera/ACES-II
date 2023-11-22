subroutine mkvector(matrix,vector,idim1)
use indices
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(in)::matrix
double precision,dimension(:),intent(inout)::vector
integer,intent(in)::idim1
! idim1 is dimension of square matrix
! idim2 is length of vector
jindex=0
do i=1,idim1
do j=1,i
jindex=jindex+1
vector(jindex)=matrix(j,i)
end do
end do

end subroutine mkvector
