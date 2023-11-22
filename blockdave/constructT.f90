subroutine constructT(newRvecs,oldRvecs,Nsize,Nblocks,&
                     totVecs,itr,maxitr,fullT,eVecs,eVals)
  integer, intent(in)::Nsize,Nblocks,itr,maxitr,totVecs
  double precision, intent(inout)::newRvecs(Nsize*totVecs)
  double precision, intent(inout)::oldRvecs(Nsize*totVecs)
  double precision,intent(inout)::fullT(totVecs,totVecs)
  double precision,intent(inout)::eVecs(totVecs,totVecs)
  double precision,intent(inout)::eVals(totVecs)

  double precision::wr(totVecs)
  double precision:: wi(totVecs),vl(1,totVecs)
  double precision:: vr(totVecs,totVecs)
  double precision::work(4*totVecs)
  integer::i,info,aoff,aend,boff,ioff,const
  double precision::total
! **Variables:
! newRvecs - R vectors fresh out of hbarxc
! Nsize    - row dim. of newRvecs
! Nblocks  - col. dim. of newRvecs
! itr      - current iteration of the block algo
! maxitr   - maximum number of iterations in block algo
! fullT    - full (albeit mostly empty) space reserved for depositing relevant
!            elements of T. 
! ScrAlpha -  temporary scratch space to store alpha
! ScrBeta  - temporary scratch space to store beta
! eVecs    - Eigenvectors of matrix T
! eVals    - Eigenvalues of matrix T
!
! **Purpose:
! Performs the block vector-block vector multiplication to create the leements
! of T, namely alpha_i=Rnew^t(Hbar*Rnew),beta_i+1=Rold^t(Hbar*Rnew)
! beta_i+1^t=Rnew^t(Hbar*Rold). Each block is of size (Nblocks x Nblocks)
! Uses current itr to determine placement; must create tridiagonal T. 
! Then Diagonalizes this mini matrix and returns eVecs/values.

  ScrAlpha=0.0d0
  ScrBeta=0.0d0
  ioff=1
  print*,'inside constructT',totVecs
  eVals=0.0d0
  eVecs=0.0d0
  do i=1,totVecs
    call getlst(newRvecs(ioff),i,1,1,1,498)
    call getlst(oldRvecs(ioff),i,1,1,1,497)
    ioff=ioff+Nsize
  enddo
  fullT=0.0d0
  const=totVecs!Nblocks*itr
  print*,'const:',const
  call dgemm('T','N',const,const,Nsize,1.0d0,oldRvecs,Nsize,&
                     newRvecs,Nsize,0.0d0,fullT,const)


#ifdef _DEBUG_LVLM
  print*,'total vecs:',totVecs
  print*,'maxiter:',maxitr
  print*,'Nblocks:',Nblocks
  print*,'size of eVals:',size(eVals)
  print*,'size of fullT:',size(fullT)
  print*,'shape of fullT:',shape(fullT)
     !call checksum("fullT stored check:",fullT,(const*const),s)
     !print*,'T matrix:'
!     call output(fullT,1,const,1,const,const,const,1)
!     call output(ScrAlpha,1,Nblocks,1,Nblocks,Nblocks,Nblocks,1)
#endif


     call dgeev("N","V",const,fullT,const,eVals,wi,vl,1,eVecs,const,work,&
                4*const,info)


#ifdef _DEBUG_LVLM
     print*,'shape of Tvecs:',shape(eVecs),size(eVecs)
     !do i=1,const
     !  call checksum("Tvecsb4sort:",eVecs(:,i),Nblocks,s)
     !enddo
#endif


     if (info.ne.0) then
        print*,'Something went wrong diagonalizing T, in constructT.f90'
        stop
     endif


#ifdef _DEBUG_LVLM
     print*,'eigenvals in T:', eVals!wr
     !call checksum("eVecs b4 sort:",eVecs,Nblocks*Nblocks,s)
#endif

     call sort(Evals,eVecs,totVecs)


#ifdef _DEBUG_LVLM
     do i=1,const
       call checksum("TvecsAFTERsort:",eVecs(:,i),totVecs,s)
     enddo
     print*,'sorted evals in T:', eVals
     print*, 'sorted Tvecs:'
!     call output(eVecs,1,const,1,const,const,const,1)
     total=0.0d0
     do i=1,const
        total=total+fullT(i,i)
     enddo
     print*,'trace of fullT:', total
     print*,'sum of eigenvalues:',sum(Evals)
#endif

end subroutine
