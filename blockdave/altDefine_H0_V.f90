subroutine altDefine_H0_V(CISmat,NSize,projP,P,initEval,H0,V)
  integer,intent(in)::NSize
  double precision,intent(in)::CISmat(NSize,NSize)
  double precision,intent(in)::projP(NSize,NSize),P(NSize,1)
  
  double precision,intent(inout)::H0(NSize,NSize),V(NSize,NSize)
  double precision,intent(inout):: initEval(1,1)
  integer::i
  initEval=matmul(matmul(transpose(P(:,1:1)),CISmat),P)
  print*,'initial guess to lowest Eval:',initEval




  H0=0.0d0
  V=0.0d0
  H0=matmul(matmul(projP,CISmat),projP)
  V=CISMAT-H0

end subroutine
