subroutine defineD(D,initEval,NSize,NBlock,H0,q,RSPT,BWPT)
  integer,intent(in)::NSize,NBlock
  double precision,intent(in)::initEval(Nblock,Nblock),H0(NSize,NSize)
  double precision,intent(in)::q(NSize,NSize-NBlock)
  double precision,intent(inout)::D(NSize-NBlock,NSize-NBlock)
  logical, intent(in):: RSPT,BWPT

  double precision::QH0Q(NSize-NBlock,NSize-NBlock)
  double precision:: tmp,tmp2
  integer::i


if (RSPT) then
  QH0Q=0.0d0
  QH0Q=matmul(transpose(q(:,:)),matmul(H0,q))
!  tmp=1.0d0/initEval(1,1)
!  print*,'tmp in D:',tmp
  D=0.0d0
  print*,'QH0Q part of Resolvent:'
  do i=1,NSize-Nblock
    print*,QH0Q(i,i)
    D(i,i)=1.0d0/(initEval(1,1)-QH0Q(i,i))
    print*,'value of D:',D(i,i)
  enddo

endif

if (BWPT) then
  QH0Q=0.0d0
  QH0Q=matmul(transpose(q(:,:)),matmul(H0,q))
  D=0.0d0
  do i=1,NSize-Nblock
    tmp2=initEval(1,1)-QH0Q(i,i)
    if (abs(tmp2).lt.1.0D-2) then
      D(i,i)=1.0d0/sign(1.0D-2,tmp2) ! changed 4.29
    else
      D(i,i)=1.0d0/(initEval(1,1)-QH0Q(i,i))
    endif
  enddo

endif


end subroutine



