       subroutine davidsonDrive(Hbar,C,NSize,cROOTStotal,NUMSOL,&
                                MAXITER,cVALS, cEVECS,TOL)
       implicit none
       integer, intent(in)::NSize,MAXITER
       integer, intent(inout)::cROOTStotal,NUMSOL
       double precision, intent(inout)::Hbar(NSize,NSize),C(NSize,NSize)
       double precision,intent(inout)::cVALS(NUMSOL)
       double precision,intent(inout)::cEVECS(Nsize,NUMSOL)
       double precision, intent(in)::TOL

       double precision::thetaOld(NUMSOL),thetaNew(NUMSOL),& 
                               diff(NUMSOL)
       double precision::s,tmp,norm,ddot
       double precision:: TVec(Nsize,Nsize)
       integer ::i,j ,I020,I030,zz
       real::time,start,finish
       logical::mask(NUMSOL)
! NSize- dimension of Hbar
! NUMSOL -Number of initial guess vectors, and increment by which we
!           expand the subspace of Hbar each loop iteration
! cROOTStotal - TOTAL NUMBER OF ROOTS CONVERGED 
! MAXITER - Maximum number of Davidson iterations--should be user defined
! Hbar - Hermitian matrix to be diagonalized
! C    - Coefficient vector meant to span a subspace of Hbar. Initialized
!        as columns of the identity matrix. Form could be improved for memory
!        considerations. In theory, only needs (Nsize x NUMSOL*MAXITER) of 
!        memory allocated, instead of (NSize x Nsize)
! Tvec - Eigenvector of size (i x i) of the tridiagonal matrix, where i is the
!        loop increment
! thetaOld-prior iterations eigenvalues
! thetaNew- current iterations eigenvalues
! cVALS   -Converged eigenvalues
! cEVECS  -Converged eigenvectors
! TOL     -Convergence tolerance
! 
! Purpose:: Block Davidson driver. Given NUMSOL columns of the identity matrix
! as guesses, iteratively expands the subspace of Hbar. Begins by performing a
! QR decomposition on the current iteration's guess space, C, and sets C=Q. Then
! uses this updated C to build a tridiagonal matrix, T, diagonalize it to find
! eigenvalues/vectors (thetaNew, TVec respectively). Using these values, builds
! a new set of guesses, r1 then q1, which are then inserted into the next rectangular block
! of C. The convergence criteria is measured by the norm of the difference b/t current and
! and prior iteration's eigenvalues. Once this falls below 0.0009, convergence
! is achieved. 
! Initialize Random Eigenvalues of size 'NUMSOL' for convergence purposes
! =>>> Set them all to 0
       I020=MAXITER*NUMSOL ! j index of C for the start of stored evals
       I030=I020+NUMSOL    ! j index of C for the start of stored evecs
       thetaOld=0.0d0
       TVec=0.0d0
       thetaNew=0.0d0
       
       time=0.0d0

       Do i=NUMSOL, MAXITER,NUMSOL
       call cpu_time(start)
       call createQ(C(:,1:i),Nsize,(i/NUMSOL)*NUMSOL)
#ifdef _DEBUG_LVLM
       call cpu_time(finish)
       print*,'Time spent is QR routine:', finish-start
       time=time+(finish-start)
! compute Matrix-vector multiplication
! T=C^t*Hbar*C; return eigenvalues/vectors to loop
       call cpu_time(start)
#endif
       call ctxhbarxc(Hbar,C(:,1:i),thetaNew,TVec(1:i,1:i),  &
                                   Nsize,i,NUMSOL)
#ifdef _DEBUG_LVLM
       call cpu_time(finish)
       print*,'Time spent is Hbarxc routine:', finish-start
       time=time+(finish-start)
! Build current iteration's Hbar eigenvecs, and extend the dimension of the search subspace

       call cpu_time(start)
#endif
       call extendCspace(Hbar,NSize,thetaNew,NUMSOL,  &
                               C(:,1:2*i),i,TVec(1:i,1:i),cEVECS)
#ifdef _DEBUG_LVLM
       call cpu_time(finish)
       print*,'Time spent is extendCspace routine:', finish-start
       time=time+(finish-start)
#endif
! Calculate difference b/t successive iterations theta to determine convergence

       diff=abs(thetaNew-thetaOld)
       print*,'List of Residual Norms for iteration ', i/NUMSOL,': '
       print*,diff
       !
!       norm=ddot(NUMSOL,diff,1,diff,1) 
!       norm=norm**(0.5d0)
!       print*,'norm:',norm
       if (all(diff.lt.TOL)) then
        print*,"** Converged on iteration: ",i/NUMSOL
        print*
        print*
        print*,'** List of converged roots: **'
        print*,thetaNew
        cVALS=thetaNew
        print*,'** in (eV) **'
        print*,thetaNew*27.2114
        cROOTStotal=NUMSOL       
        print*,'** Orthonormalizing vectors:'
        call GramSchmidt(cEVECS,NSize,NUMSOL)
        tmp=0.0d0
        do zz=1,NUMSOL
          tmp=dot_product(cEVECS(:,zz),cEVECS(:,zz))
          tmp=1.0d0/sqrt(tmp)
          cEVECS(:,zz)=tmp*cEVECS(:,zz)
          print*,'factor on denom',tmp
          call checksum("finalCISVec:",cEVECS(:,zz),NSize,s) 
        enddo
        exit
       endif
       
       thetaOld=thetaNew
       print*,'theta old',thetaOld*27.2114
       enddo 
        mask=diff.lt.TOL
        print*,'Number of roots satisfying convergence criteria:'&
                                                        ,count(mask)
       cROOTStotal=count(mask)
#ifdef _DEBUG_LVLM
       print*,'** Total time for block (Hermitian) Davidson: ',time/3600 
#endif       
       return
       end
