subroutine redPartDiag(CISmat,NSize,TOL)
  integer,intent(in)::Nsize
  double precision,intent(in)::TOL
  double precision,intent(in)::CISmat(NSize,NSize)
  double precision::P(NSize,1),Q(Nsize,NSize-1)


  double precision::projP(NSize,NSize),projQ(NSize,NSize)
  double precision::Pert(NSize,NSize),Denom(NSize,NSize)
integer::nbas,natm,nocc,nvirt,scratch(2)
  integer::order

  double precision::H0(NSize,NSize),tmpHold(NSize-1,NSize-1)
  double precision::V(NSize,NSize),D(NSize-1,NSize-1)
  double precision::V0out(NSize,NSize),pVq(1,NSize-1)
  double precision::qVq(NSize-1,NSize-1)
  double precision::E_0(1,1),E_1(1,1)
  double precision,allocatable::AA(:,:),BB(:,:),CC(:,:)
  double precision::tmp(1,1),initEval(1,1),evaltmp(1,1)
  double precision::guessE(1,1),correction(1,1)
  logical::BWPT,RSPT,Converged,sVecAlgo
  

  integer::NBlock
!      call getrec(1,'JOBARC','NBASTOT',1,nbas)
!      call getrec(1,'JOBARC','NREALATM',1,natm)
!        call getrec(1,'JOBARC','NOCCORB',2,scratch)
!      nocc=scratch(1)
!      nvirt=nbas-nocc
!print*,'nocc:',nocc
!print*,'nvirt:',nvirt
!print*,'NSize:',NSize,nocc*nvirt
!print*,'Nbas:',nbas



  BWPT=.True.
  RSPT=.False.


!if (RSPT) then
!  call defineD(D,initEval,NSize,H0,Q,guessE,RSPT,BWPT)
!
!  xVq=0.0d0
!  qVDx=0.0d0
!  qVDq=0.0d0
!! RS PT correction:
!
!  tmpHold=0.0d0
!  evaltmp=initEval(1,1)
!  do order=0,1
!    allocate(AA(order,order),BB(order,order),CC(order,1))
!
!    call buildRSPert(V,D,P,Q,NSize,pVq,qVq,tmpHold,order,E_0,evaltmp)
!    print*
!    print*,'****************************'
!    print*,'Order: ', order+1
!    print*,'Correction at this order:', E_0
!    evaltmp=evaltmp+E_0(1,1)
!    print*,'Current excitation energy:',evaltmp
!    print*,'****************************'
!
!    deallocate(AA,BB,CC)
!  enddo
!
!endif

if (BWPT) then
  NBlock=1
  print*,'Entering partitioning algorithm....'
!  call partitionAlgo(CISmat,NSize,NBlock)

  print*,'**********************************'
  print*,'**********************************'
  print*,'**********************************'
  print*,'**********************************'
  print*,'**********************************'
  print*,'**********************************'
  print*,'Entering single-vector variational algorithm now....'
  sVecAlgo=.True.
!  call variationalAlgo(CISmat,NSize,NBlock,sVecAlgo)


  print*
  print*,'Entering block-variational algorithm now...'
  NBlock=5
  sVecAlgo=.False.
!  call variationalAlgo(CISmat,NSize,NBlock,sVecAlgo)





 print*,'Entering Pade approximant algo now....'
 Nblock=1
 sVecAlgo=.True.
 call padeAlgo(CISmat,NSize,NBlock,sVecAlgo)

  !call blockVariationalAlgo(CISmat,NSize,NBlock,sVecAlgo)
!  guessE=initEval(1,1)
!  print*,'Initial guess:', guessE
!  pVq=0.0d0
!  qVq=0.0d0
!
!  qVDx=0.0d0
!  qVDq=0.0d0
!! RS PT correction:
!
!  tmpHold=0.0d0
!  evaltmp=initEval(1,1)
!  do order=0,0
!
!    call defineD(D,initEval,NSize,H0,Q,guessE,RSPT,BWPT)
!    call buildBWPert(V,D,P,Q,NSize,pVq,qVq,tmpHold,order,E_0)
!    print*
!    print*,'****************************'
!    print*,'Order: ', order+1
!    print*,'Correction at this order:', E_0
!    evaltmp=evaltmp+E_0(1,1)
!    print*,'Current excitation energy:',evaltmp
!    print*,'****************************'
!    guessE=evaltmp
!    
!
!
!
!  enddo

endif


end subroutine  
