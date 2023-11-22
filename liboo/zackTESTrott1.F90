      subroutine zackTESTrott1(nocc,nvirt)
      IMPLICIT INTEGER (A-Z)
      integer,intent(in)::nocc,nvirt
      DOUBLE PRECISION::iden(nocc+nvirt,nocc+nvirt)
      double precision,allocatable::scfVecs(:,:),kappaVO(:,:)
      double precision,allocatable::kappa(:,:)
      double precision,allocatable::rotatedVecs(:,:)
      double precision, allocatable::newSCFvirts(:,:)
      double precision, allocatable::newSCFocc(:,:)
      double precision, allocatable::AOovrlap(:,:)
      double precision, allocatable::AOtrans(:,:)
      double precision, allocatable::scr(:)
      double precision, allocatable::AOtranEvecs(:,:)
      double precision, allocatable::finalSCFvec(:,:)
      double precision,allocatable::transVec(:,:)
      double precision, allocatable::Xtrans(:,:)
      double precision, allocatable::fullGrad(:,:)
        integer ::max_step
      integer::nbas,vecDim,Naobfns,nbas2
        double precision::sub_grad(nvirt,nocc)
      double precision::scftmps(nocc+nvirt,nocc+nvirt)
      logical::converge

      print*,'@zackRott1: insides zack rot'
      nbas=nocc+nvirt
      nbas2=nbas**2

      allocate(scfVecs(nbas,nbas),kappaVO(nvirt,nocc))
      allocate(kappa(nbas,nbas),rotatedVecs(nbas,nbas))      
      allocate(newSCFvirts(nbas,nbas))
      allocate(newSCFocc(nbas,nbas),AOovrlap(nbas,nbas))
      allocate(AOtrans(nbas,nbas),scr(2*nbas2))
      allocate(AOtranEvecs(nbas,nbas))
      allocate(finalSCFvec(nbas,nbas))
      allocate(transVec(nbas,nbas),Xtrans(nbas,nbas))
      allocate(fullGrad(nbas,nbas))
      CALL GETREC(20,'JOBARC','SCFEVECA',NBAS2,scfVecs)
      call checksum('@zacRott1: checksum postscfVecs:',scfVecs,nbas2)
      call getrec(20,'JOBARC','SCFEVCA0',NBAS2,scftmps)
      call checksum('@zacRott1: checksum prescfVecs:',scftmps,nbas2)
      ! Load orbital gradient, check convergence
      call getrec(20,'JOBARC',"ORBRTGRD",nbas2,fullGrad)
      print*,'@zackTESTrott1: checking RMSE of gradient'
      call checkGrad(fullGrad,nbas2,10E-5,converge)

      ! Extract only the vo part of full orbital gradient
      ! and rescale each element by diagonal elements of Fock matrix
      sub_grad=fullGrad(nocc+1:nbas,1:nocc)
      call ScaleKappaOV(kappaVO,sub_grad,nocc,nvirt,max_step)

      ! Stored rescaled OV part of kappa in the appropriate blocks
      ! of the full kappa
      kappa=0.0d0
      kappa(1:nocc,nocc+1:nbas)=-1.0d0*transpose(kappaVO)
      kappa(nocc+1:nbas,1:nocc)=kappaVO


       ! **NOTE** uncomment next 2 lines to include oo, vv 
       ! (unscaled) parts of gradient into full Kappa
!      kappa(1:nocc,1:nocc)=-1.0d0*fullGrad(1:nocc,1:nocc)
!      kappa(nocc+1:nbas,nocc+1:nbas)=-1.0d0*fullGrad(nocc+1:nbas,&
!     &                                          nocc+1:nbas)


      ! Now, determine the new set of MO coefficients using 
      ! old SCF MO C and the rotation matrix, kappa
      ! 
      rotatedVecs=0.0d0
      rotatedVecs=matmul(scfVecs,kappa) + scfVecs
      call checksum('@zacRott1: checksum rotatedVec',rotatedVecs,nbas2)
      call putrec(20,'JOBARC','SCFEVECA',nbas2,rotatedVecs)



      ! load AO ovrlap matrix and transform to MO basis
      call getrec(20,'JOBARC','AOOVRLAP',nbas2,AOovrlap)
      call AO2MO2(AOovrlap,AOovrlap,rotatedVecs,scr,nbas,nbas,1)
      CALL CHECKSUM('@zacRott1: checksum trans.Covrlap',AOovrlap,nbas2)
!      print*,'output of trans. overlap'
!      call output(AOovrlap,1,nbas,1,nbas,nbas,nbas,1)


      ! Calculate X=U s^(-1/2)U^t via lowdin orthogonalization
      call eig(AOovrlap,rotatedVecs,nbas,nbas,0)
      call invsqt(AOovrlap,nbas+1,nbas)
      call checksum('@zacRott1:chksum Smat eig vecs:',rotatedVecs,nbas2)
      call checksum('@zacRott1:chksum Smat eig vals:',AOovrlap,nbas2)
      Xtrans=0.0d0
      Xtrans=matmul(matmul(rotatedVecs,AOovrlap),transpose(rotatedVecs))
      print*,'canonical ortho. trans.'
      call checksum('@zacRott1: chksum otho.trans.',Xtrans,nbas2)


      ! NOW TRANSFORM MO VECTOR TO ORTHOGONAL REPRESENTATION AND DEPOSIT
      ! C =C' * X, where S has been trans. to MO basis
      ! Alternatively, we could have done
      call getrec(20,'JOBARC','SCFEVECA',NBAS2,finalSCFvec)
      transVec=0.0d0
      transVec=matmul(finalSCFvec,Xtrans)
      call checksum('@zacRott1:chksum transVec putrec',transVec,nbas2)
      call putrec(20,'JOBARC','SCFEVCA0',nbas2,transVec)



      ! Check that C'^T S C' ==1
        call getrec(20,'JOBARC','AOOVRLAP',nbas2,AOovrlap)
        iden=matmul(matmul(transpose(transVec),AOovrlap),transVec)
        print*,'@zacRott1: Checking that C^T S C ==1'
        call output(iden,1,nbas,1,nbas,nbas,nbas,1)

#ifdef _NO_SKIP
      print*,'Eigenvalues of S^-1/2:'
      call output(AOtrans,1,nbas,1,nbas,nbas,nbas,1)
      print*,'Eigenvectors of S^-1/2:'
      call output(rotatedVecs,1,nbas,1,nbas,nbas,nbas,1)
#endif 
      deallocate(scfVecs,kappaVO,kappa,rotatedVecs,newSCFvirts) 
      deallocate(AOtranEvecs,scr,newSCFocc,AOovrlap,AOtrans)
      deallocate(fullGrad,finalSCFvec,transVec,Xtrans)




      end subroutine
