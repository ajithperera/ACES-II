subroutine defineH0(CISmat,NSize,NBlock,P,initEval,H0)
  integer,intent(in)::NSize,Nblock
  double precision,intent(in)::CISmat(NSize,NSize)
  double precision,intent(in)::P(NSize,NBlock)
  
  double precision,intent(inout)::H0(NSize,NSize)
  double precision,intent(inout):: initEval(NBlock,NBlock)
  integer::i
  print*,'p'
  print*,p
!  initEval=matmul(matmul(transpose(P(:,1:NBlock)),CISmat),P)
!  print*,'initial guess to lowest Eval:',initEval




  H0=0.0d0
  print*,'Diag of H0'
  do i=1,NSize
    H0(i,i)=CISmat(i,i)
    print*,H0(i,i)
  enddo
  initEval=matmul(matmul(transpose(P(:,1:NBlock)),H0),P)
  print*,'initial guess to lowest Eval:',initEval
print*,'H0'
call output(H0,1,NSize,1,NSize,NSize,NSize,1)
end subroutine
