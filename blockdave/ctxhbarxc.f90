      subroutine ctxhbarxc(Hbar,C,thetaNew,eVec,NSize,i,initVec)
      integer,intent(in)::NSize,i,initVec
      double precision, intent(in):: Hbar(Nsize,Nsize),C(Nsize,i)
      double precision, intent(inout)::thetaNew(initVec)
      double precision, intent(inout)::eVec(i,i)
      double precision:: T(i,i),Scr(Nsize,i)
      integer::a,junk,info
      double precision::wr(i),wi(i),vl(1,i),vr(i,i),work(4*i)
      integer::resortedArr(i)
      logical,dimension(i)::mk
      double precision:: thetaCPY(i),eVecCPY(i,i)
! Purpose: Computes the matrix-vector multiplication of C^t*Hbar*C
!          which is stored in T. T is then diagonalized to form evecs and
!          evalues.
! In general, initVec < i unless we are in first iteration of Davidson loop

      T=0.0d0

! vector-Matrix-vector multiplication

      call dgemm('N','N',NSize,i,Nsize,1.0d0,Hbar,Nsize,C,Nsize,0.0d0,Scr,Nsize)

!call dgemm('T','N',i,i,Nsize,1.0d0,C,i,Scr,i,0.0d0,T,i)

      T=matmul(transpose(C),Scr)

! Diagonalize tridiagonal matrix T, then copy evalues to thetaNew
! This routine works because T is symmetric; otherwise use blas then sort
! the resulting evecs/values

      !call eig(T,vr,junk,i,0)
!      do a=1,initVec
!         wr(a)=T(a,a)
!      enddo
      print*,'Using new dgeev'
      call dgeev('N','V',i,T,i,wr,wi,vl,1,vr,i,work,4*i,info)
      resortedArr=0
      mk=.true.
      do k=1,i
        resortedArr(k)=minloc(wr,1,mk)
        mk(minloc(wr,mk))=.false.
      enddo
      thetaCPY=0.0d0
      eVecCPY=0.0d0
      do k=1,i
        place=resortedArr(k)
        thetaCPY(k)=wr(place)
        eVecCPY(:,k)=vr(:,place)
      enddo
        thetaNew=thetaCPY(:initVec)
        eVec=eVecCPY
        
!      do a=1,initVec
!         thetaNew(a)=T(a,a)
!      enddo
#ifdef _DEBUG_LVLM
      print*,'thetaNew in ctxhbarxc'
      print*,thetaNew
#endif
!  print*,'evals:'
!  print*,T
!  print*,'evecs'
!  print*,eVec
!  print*,'matrix mult evec^t T evec'
!  print*,matmul(transpose(eVec),matmul(T,eVec))

   Return
   end 
