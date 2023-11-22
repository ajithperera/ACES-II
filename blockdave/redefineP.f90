subroutine redefineP(p,fullSpace,R0,V,pNew,NSize,NBlock,itr)
  integer,intent(in)::NSize,Nblock,itr
    double precision,intent(in)::R0(NSize,NSize),V(NSize,NSize)
  double precision,intent(in)::p(NSize,NBlock),fullSpace(NSize,Nblock*100)
  double precision,intent(inout)::pNew(NSize,NBlock)
  integer::i
  double precision::tmp,bothPs(NSize,Nblock*2)
  bothPs=0.0d0
!  pNew=0.0d0
  pNew(:,1)=p(:,1)+fullSpace(:,1+itr)

! Normalize pNew
  do i=1,NBlock
    tmp=dot_product(pNew(:,i),pNew(:,i))
    tmp=1.0d0/sqrt(tmp)
    print*,'Normalizing weight for pNew:',tmp
 !   pNew(:,i)=pNew(:,i)*tmp
  enddo

! Do I orthogonalize pNew against p...?
 print*,'Overlap b/t p and pNew:',dot_product(p(:,1),pNew(:,1))

! bothPs(:,1:Nblock)=p(:,:NBlock)
! bothPs(:,NBlock+1:2*NBlock)=pNew(:,:NBlock)
! call GramSchmidt(bothPs,NSize,NBlock*2) 

! orthogonalized pNew is on 2nd col. of bothPs
! pNew=bothPs(:,Nblock+1:2*NBlock)
! print*,'Overlap b/t p and pNew after GS:',dot_product(p(:,1),pNew(:,1))


end subroutine
