subroutine defineV(CISmat,NSize,V)
  integer,intent(in)::NSize
  double precision,intent(in)::CISmat(NSize,NSize)
  double precision,intent(inout)::V(NSize,NSize)

  integer::i,MATDIM
  V=0.0d0 
  MATDIM=NSize
  V=CISmat
  do i=1,NSize
    V(i,i)=0.0d0
  enddo
print*,'V'
call output(V,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
end subroutine
