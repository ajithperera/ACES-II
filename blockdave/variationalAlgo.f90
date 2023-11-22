subroutine variationalAlgo(CISmat,NSize,Nblock,sVecAlgo)
  integer,intent(in)::NSize,Nblock
  double precision,intent(in)::CISmat(NSize,NSize)
  logical,intent(in)::sVecAlgo

  double precision::Hbar(NSize,NSize),P(NSize,NBlock),Q(Nsize,NSize-NBlock)
  double precision::projP(NSize,NSize),projQ(NSize,NSize)
  logical::Converged,BWPT,RSPT
  integer::offset,i,j,m,z,iter,order


  double precision::H0(NSize,NSize),V(NSize,NSize)
  double precision::R0(NSize*NSize*NBlock)
  double precision:: tse(NBlock,NBlock),E0(NBlock,NBlock)
  double precision::eps0(NBlock,NBlock),tseOld(NBlock,NBlock)
  double precision::pNew(NSize,NBlock)
  double precision::fullSpace(NSize,NBlock*100)

  tseNew=0.0d0
  pNew=0.0d0
  fullSpace=0.0d0
  tseOld=0.0d0
  Hbar=CISmat
  BWPT=.TRUE.
  RSPT=.FALSE.
  Converged=.False.
  call definePQ(CISmat,NSize,Nblock,P,Q,projP,projQ,sVecAlgo)

  iter=1
  order=0
  do while (.not.Converged)
    print*
    print*,'************************'
    print*,'************************'
    print*,'Beginning iteration:',iter
    if (order.eq.0) then
      call defineH0(Hbar,NSize,Nblock,P,tse,H0)
      call defineV(Hbar,NSize,V)
      print*,'Guess evalue:',tse
    endif

    if (sVecAlgo) then ! If one root at a time
      z=1
      call defineR0(R0,tse,NSize,Nblock,H0,q,RSPT,BWPT,sVecAlgo,z)
      call extPSpace(p,R0,V,NSize,NBlock,order,iter,fullSpace,sVecAlgo,z)
    else               ! If block Algo is used
      offset=1
      do m=1,NBlock
        call defineR0(R0(offset),tse,NSize,Nblock,H0,q,RSPT,BWPT,sVecAlgo,m)
        call extPSpace(p,R0(offset),V,NSize,NBlock,order,iter,fullSpace,sVecAlgo,m)
        offset=offset+NSize*NSize
      enddo
    endif

    call buildRedH(Hbar,tse,fullSpace,NSize,Nblock,order,iter,sVecAlgo)
    print*,'Current root:', tse
    if (abs(abs(tse(1,1))-abs(tseOld(1,1))).gt.1.0D-10) then ! not converged
      tseOld=tse 
    else ! all roots are converged
      print*,'Converged on iteration:',iter
      print*,'Converged at BWPT order:',order
      print*,'Converged root:', tse
      Converged=.True.
    endif

    order=order+1
    iter=iter+1

  enddo

end subroutine
