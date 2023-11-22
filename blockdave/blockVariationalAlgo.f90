subroutine blockVariationalAlgo(CISmat,NSize,Nblock)
  integer,intent(in)::NSize,Nblock
  double precision,intent(in)::CISmat(NSize,NSize)


  double precision::Hbar(NSize,NSize),P(NSize,NBlock)
  double precision::Q(Nsize,NSize-NBlock)
  double precision::projP(NSize,NSize),projQ(NSize,NSize)
  logical::sVecAlgo,Converged,BWPT,RSPT
  integer::i,j,iter,order


  double precision::H0(NSize,NSize),V(NSize,NSize)
  double precision::R0(NSize,NSize)
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
  sVecAlgo=.False.
  call definePQ(CISmat,NSize,1,P,Q,projP,projQ,sVecAlgo)


end subroutine
