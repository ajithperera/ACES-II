subroutine defineR0(R0,tse,NSize,Nblock,H0,q,RSPT,BWPT,sVecAlgo,m)
  integer,intent(in)::NSize,Nblock,m
  double precision,intent(in)::tse(NBlock,NBlock),H0(NSize,NSize)
  double precision,intent(in)::q(NSize,NSize-Nblock)
  double precision,intent(inout)::R0(NSize,NSize)
  logical, intent(inout):: RSPT,BWPT,sVecAlgo

  double precision::Denom(NSize-Nblock,NSize-Nblock)
  integer::MATDIM
  double precision::tseTmp(Nblock,Nblock)
MATDIM=NSize-Nblock
Denom=0.0d0
tseTmp=0.0d0
if (BWPT) then
  if (sVecAlgo) then
    call defineD(Denom,tse,NSize,Nblock,H0,q,RSPT,BWPT)
    R0=matmul(matmul(q,Denom),transpose(q(:,:)))
    !print*,R0
    print*,'printing Denom'
    call output(Denom,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
    print*,'printing negative definite R0'
    call output(R0,1,NSize,1,NSize,NSize,NSize,1)
   else
     tseTmp(1,1)=tse(m,m)
     print*,'Using tse:',tseTmp(1,1),'on R0^i',m
     call defineD(Denom,tseTmp,NSize,Nblock,H0,q,RSPT,BWPT)
     R0=matmul(matmul(q,Denom),transpose(q(:,:)))
     print*,'printing Denom'
     call output(Denom,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
     print*,'printing negative definite R0'
     call output(R0,1,NSize,1,NSize,NSize,NSize,1)
   endif
else
  print*,'Havent written code for RSPT yet....'
endif
end subroutine
