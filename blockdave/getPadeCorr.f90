subroutine getPadeCorr(CISmat,NSize,Nblock,sVecAlgo,corr,order,R0,p)
  integer,intent(in)::NSize,Nblock,order
  double precision,intent(in)::CISmat(NSize,NSize)
  logical,intent(in)::sVecAlgo
  double precision,intent(inout)::corr(Nblock,Nblock)
  double precision,intent(in)::R0(Nsize,NSize),p(NSize,Nblock)

  double precision::matA(order*Nblock,order*Nblock)
  double precision::matB(order*Nblock,order*Nblock)
  double precision::vecC(order*Nblock,Nblock)
  double precision::tempEps(Nblock,Nblock)
  double precision::combAB(order*Nblock,order*Nblock)
  double precision:: combABbkup(order*Nblock,order*Nblock)
  integer::i,j,n,ioff,joff,offset,info
  double precision:: work(Nblock*order)
  integer::ipiv(NBlock*order)


  combAB=0.0d0
  combABbkup=0.0d0
! Performs C^t (A-B)^(-1) C to obtain correction to initial guess
! To do this, builds vecC, matA, matB which depend on the order. 

! Build A:
  ioff=1
  joff=1
  offset=Nblock-1
  do i=0,order-1
    do j=0,order-1
      call buildEpsBlock(CISmat,NSize,Nblock,sVecAlgo,i+j,R0,p,tempEps)
    
      matA(ioff:ioff+offset,joff:joff+offset)=tempEps(:,:)
      joff=joff+Nblock
     enddo
     ioff=ioff+Nblock
     joff=1
   enddo    


! Build vecC, which is just the first column of matA:

  vecC=matA(:,:Nblock)


! Build B:
  ioff=1
  joff=1
  do i=0,order-1
    do j=1,order
      call buildEpsBlock(CISmat,NSize,Nblock,sVecAlgo,i+j,R0,p,tempEps)

      matB(ioff:ioff+offset,joff:joff+offset)=tempEps(:,:)
      joff=joff+Nblock
     enddo
     ioff=ioff+Nblock
     joff=1
   enddo


! Compute (A-B)^-1
! source:::  http://fortranwiki.org/fortran/show/Matrix+inversion
  combAB=matA-matB
  combABbkup=combAB
  n=size(combAB,1)
  call dgetrf(n,n,combABbkup,n,ipiv,info)
  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  call dgetri(n,combABbkup,n,ipiv,work,n,info)
  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if


! Compute C^t (A-B)^-1 C

  corr=matmul(transpose(vecC(:,:)),matmul(combABbkup,vecC))
  print*,'correction to initial guess is: ', corr




end subroutine
