subroutine diagFULL(Mat,dims,Nblocks,eVecs,eVals)
  integer,intent(in)::dims,Nblocks
  double precision,intent(in)::Mat(dims,dims)
  double precision, intent(inout)::eVecs(dims,dims)
  double precision, intent(inout)::eVals(Nblocks)

  double precision::wr(dims),wi(dims),vl(1,dims)
  double precision:: vr(dims,dims),work(4*dims)
  integer::info
!! **Purpose::
! Diagonalize the mini matrix T, and returns the relevant eigenvalues/vectors
! sorted in sequential order.


     call dgeev("N","V",dims,Mat,dims,wr,wi,vl,1,vr,dims,work,&
                4*dims,info)
     if (info.ne.0) then
        print*,'Something went wrong diagonalizing T, in constructT.f90'
        stop
     endif
     call sort(wr,vr,dims)

    eVecs=vr
    eVals=wr(1:Nblocks)


end subroutine
