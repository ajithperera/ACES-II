subroutine extPSpace(p,R0,V,NSize,NBlock,order,itr,fullSpace,sVecAlgo,m)
  integer,intent(in)::NSize,Nblock,order,itr,m
  double precision,intent(in)::p(NSize,Nblock)
  double precision,intent(in)::V(NSize,NSize),R0(NSize,NSize)
  double precision,intent(inout)::fullSpace(NSize,100*Nblock)
  logical,intent(in)::sVecAlgo

  double precision::R0V(NSize,NSize),tmp(Nsize,NSize)
  double precision::tmpP(NSize,1)!NBlock)
  integer::k

  R0V=matmul(R0,V)
  tmp=R0V
  if (order.eq.0 .and. itr.eq.1) then ! change here 4.29.21
    fullSpace(:,1:NBlock)=p
  endif
! build correction (R0V)^(order+1)|p>

  do k=1,order!0,order ! change here 4.29.21
    tmp=matmul(tmp,R0V)
  enddo

  if (sVecAlgo) then
    print*,'original p:',p
    print*,'R0V'
call output(R0V,1,NSize,1,NSize,NSize,NSize,1)
    tmpP=matmul(tmp,p)
    print*,'correction to p:',tmpP
!  Add correction to full space |R>
!  if (sVecAlgo) then ! **SINGLE VECTOR ONLY
    fullSpace(:,(itr)*Nblock+1:(1+itr)*Nblock)=tmpP

! Orthogonalize space: |R>=|p phi0 phi1 ....phip>
! See eq. 
    call GramSchmidt(fullSpace,NSize,(1+itr)*Nblock)
    print*,'orthogonalize correction to p:',fullSpace(:,:(1+itr)*Nblock)


  else !! block algo---takes one R0^(m) at a time  
    tmpP(:,1)=matmul(tmp,p(:,m))
    fullSpace(:,itr*Nblock+m)=tmpP(:,1)
    call GramSchmidt(fullSpace,NSize,itr*Nblock+m)
!    tmpP=matmul(tmp,p) !(:,m))
!    fullSpace(:,(itr)*Nblock+1:(1+itr)*Nblock)=tmpP
!    call GramSchmidt(fullSpace,NSize,(1+itr)*Nblock)

  endif


end subroutine
