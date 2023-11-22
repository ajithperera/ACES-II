subroutine buildBWPert(V,D,p,q,Nsize,pVq,qVq,tmpHold,k,E_k)
  integer,intent(in)::k,Nsize
  double precision,intent(in)::p(NSize,1),q(NSize,NSize-1)
  double precision,intent(in)::V(NSize,NSize),D(NSize-1,NSize-1)
  double precision,intent(inout)::pVq(1,NSize-1),qVq(NSize-1,NSize-1)
  double precision,intent(inout)::E_k(1,1)
  double precision,intent(inout)::tmpHold(NSize-1,NSize-1)


!  double precision::tmpHold(NSize-1,Nsize-1)


  double precision:: tmpV(NSize,NSize)

  integer::i,j,MATDIM



  MATDIM=NSize-1
  if (k.eq.0) then
    pVq=matmul(transpose(p(:,1:1)),matmul(V,q)) 
    E_k=matmul(matmul(pVq,D),transpose(pVq(:,:)))

  else if (k.eq.1) then
    qVq=matmul(transpose(q(:,:)),matmul(V,q))
    tmpHold=matmul(D,matmul(qVq,D))
    print*,'tmpHOld:'
    call output(tmpHold,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
    print*
    print*,'pVq'
    call output(pVq,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
    E_k=matmul(matmul(pVq,tmpHold),transpose(pVq(:,:)))


!    qVDq=matmul(transpose(q(:,:)),matmul(tmpVD,q(:,:)))
!    E_k=matmul(transpose(qVDx(:,:)),matmul(qVDq,transpose(pVq(:,:))))

  else
     tmpHold=matmul(tmpHold,matmul(qVq,D))
!     print*,'tmphold for order:', k+1
!    call output(tmpHold,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
     do i=2,k
       tmpHold=matmul(matmul(tmpHold,qVq),D)
     enddo

     E_k=matmul(matmul(pVq,tmpHold),transpose(pVq(:,:)))


!     E_k=matmul(pVq,matmul(tmpHold,transpose(pVq(:,:))))
!    tmpHold=qVDq
!    do i=1,k-1
!      tmpHold=matmul(qVDq,tmpHold)
!    enddo
!    E_k=matmul(transpose(qVDx(:,:)),matmul(tmpHold,transpose(pVq(:,:))))
  endif


    
end subroutine 
