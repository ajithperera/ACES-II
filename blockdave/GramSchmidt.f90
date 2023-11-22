subroutine GramSchmidt(A,nrows,ncols)
  integer,intent(in)::ncols,nrows
  double precision,intent(inout)::A(nrows,ncols)

  integer::j,k
  double precision::temp,norm


  !print*,'rows,cols:',nrows,ncols
  !print*,'A:',A
!print*
!print*
  do j=1,ncols
    do k=1,j-1
      !print*,'ak,aj',A(:,k),1,A(:,j)
      temp=dot_product(A(:,k),A(:,j))
     ! print*,'temp',temp
      A(:,j)=A(:,j)-temp*A(:,k)
    enddo
    !print*,'A(:,j)',A(:,j)
    norm=dot_product(A(:,j),A(:,j))
    norm=sqrt(norm)
    !print*,'norm:',norm
    A(:,j)=A(:,j)/norm
  enddo

end subroutine
    
