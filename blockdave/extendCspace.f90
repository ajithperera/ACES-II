      subroutine extendCspace(Hbar,NSize,thetaNew,initVec,C,i,Tvec,&
                                cEVECS)
      integer, intent(in)::Nsize,initVec,i
      double precision, intent(inout)::Hbar(Nsize,Nsize)
      double precision, intent(inout)::thetaNew(initVec)
      double precision, intent(inout)::Tvec(i,i)
      double precision,intent(inout)::C(NSize,2*i)
      double precision, intent(inout)::cEVECS(Nsize,initVec)
      double precision:: eye(Nsize,Nsize)
      double precision:: Scr(NSize),HbarScr(Nsize,Nsize)
      double precision:: q1(Nsize),r1(Nsize),factor
      integer::j,k

! Purpose:: Use C(:,i),Tvec with Hbar,thetaNew to expand the guess space to C(:,2*i)
! Actions:: Create Identity matrix, weight it by theta(j), subtract this
!           from Hbar, (NSize x NSize). Multiply this by the multplication of
!           C(:,i)*Tvec(i,i) ======> call result r1 of dimension (NSize x 1).
!           Normalize r1 w.r.t (thetaNew(j)-Hbar(j,j)), then store in column
!           C(:,i+j).
! Hbar - full, original matrix
! NSize - dimension of Hbar
! thetaNew - list of dimension (initVec) containing the lowest
!            eigenvalues from diagonalization of tridiagonal matrix T
! initVec - User specified number of initial guess vectors
! C - Eigenvector matrix that spans a subspace of Hbar
! i - Current iteration of the Block Davidson algo
! Tvec - Eigenvector matrix of the tridiagonal T matrix, of dimension (i x i)

      HbarScr=0.0d0

      do j=1,initVec
         Scr=0.0d0
         HbarScr=Hbar
         r1=0.0d0

! Compute Scr=C(:,1:i)*Tvec(:,j)
! In other words, compute the current iteration's Hbar eigenvectors and store
! in Scr.
!   call dgemm('N','N',NSize,1,i,1.0d0,C(:,1:i),NSize,Tvec(:,j),i,Scr,NSize)

         Scr=matmul(C(:,1:i),Tvec(:,j))
         cEVECS(:,j)=Scr

! Compute HbarScr=Hbar - thetaNew(j)*I. HbarScr has dimension (NSize x
! Nsize)
  
          eye=0.0d0
          do k=1,Nsize
            eye(k,k)=1.0d0
          enddo
          call daxpy(Nsize**2,-thetaNew(j),eye,1,HbarScr,1)


! Compute r1=(Hbar-thetaNew(j)*I) * Scr. 'r1' has dimension (Nsize x 1)

          call dgemm('N','N',NSize,1,NSize,1.0d0,HbarScr,NSize,Scr,NSize,0.0d0,r1,NSize)

! Normalize r1 to form a q1 =>>> store into the (i+j) column of C
! q1=r1/(thetaNew(j) -Hbar(j,j)) then C(:,i+j)=r1

         q1=0.0d0
         factor=1.0d0/(thetaNew(j)-Hbar(j,j))
         call daxpy(NSize,factor,r1,1,q1,1)
         C(:,i+j)=q1

      enddo

      Return
      End 

