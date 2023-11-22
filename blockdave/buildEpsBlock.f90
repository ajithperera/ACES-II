subroutine buildEpsBlock(CISmat,NSize,Nblock,sVecAlgo,k,R0,p,tempEps)
  integer,intent(in)::NSize,Nblock,k
  double precision,intent(in)::CISmat(NSize,NSize)
  logical,intent(in)::sVecAlgo
  double precision,intent(in)::R0(Nsize,NSize),p(NSize,Nblock)
  double precision,intent(inout)::tempEps(Nblock,Nblock)


  double precision::HR0H(Nsize,Nsize),R0H(NSize,Nsize)
  integer::i
! builds eps_k = <p|H(R0H)^(k+1)|p>
  HR0H=matmul(CISmat,matmul(R0,CISmat))
  R0H=matmul(R0,CISmat)
  do i=1,k
    HR0H=matmul(HR0H,R0H)
  enddo

  tempEps=matmul(transpose(p(:,:)),matmul(HR0H,p))






end subroutine
