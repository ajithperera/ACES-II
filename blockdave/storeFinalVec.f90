subroutine storeFinalVec(oldRvecs,Tvecs,NSize,Nblocks,itr,totVecs)
  integer, intent(in)::NSize,Nblocks,itr,totVecs
  double precision,intent(inout)::oldRvecs(NSize,totVecs)
  double precision,intent(in)::Tvecs(totVecs,totVecs)

  double precision:: finalVec(NSize,totVecs),overlap
  integer::i

  finalVec=matmul(oldRvecs,Tvecs)


  do i=1,Nblocks
    overlap=0.0d0
    do j=1,Nsize
      overlap=overlap+finalVec(j,i)*finalVec(j,i)
    enddo
#ifdef _DEBUG_LVLM
    print*,'finalVec:',i
    print*,'overlap:',overlap
#endif
  enddo

  do i=1,Nblocks
    oldRvecs(:,i)=finalVec(:,i)
    call putlst(finalVec(:,i),i,1,1,1,497)
  enddo





end subroutine
