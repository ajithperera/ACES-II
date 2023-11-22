subroutine contractVR0V(V,Denom,Nsize,nocc,nvirt)
 integer,intent(in)::nsize,nocc,nvirt
 double precision,intent(in)::Denom(NSize,NSize)
 double precision,intent(inout)::V(NSize,NSize)

V=matmul(V,Denom)

end subroutine
