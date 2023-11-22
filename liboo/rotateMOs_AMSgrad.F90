       subroutine rotateMOs_AMSgrad(kappa,gradVecDim,Naobfns,nocc,nvirt)
        implicit none
        integer, intent(in)::gradVecDim,Naobfns,nocc,nvirt
        double precision,intent(in)::kappa(gradVecDim,gradVecDim)

        double precision::Umat(gradVecDim,gradVecDim)
        double precision::oldMO(gradVecDim**2),newMO(gradVecDim**2)
        double precision::secondOrderX(gradVecDim,gradVecDim)
        integer::i,Isqlen,nbas,nbas2
        double precision:: checkUnitary(gradVecDim,gradVecDim)
        double precision:: tmpX(gradVecDim,gradVecDim)

        double precision:: scfVecs(gradVecDim,gradVecDim)
        double precision:: rotatedVecs(gradVecDim,gradVecDim)
        double precision::AOovrlap(gradVecDim,gradVecDim)
        double precision::scr(2*gradVecDim**2)
        double precision::xx(gradVecDim,gradVecDim)
        double precision::finalSCFvec(gradVecDim,gradVecDim)
        double precision::transVec(gradVecDim,gradVecDim)
        double precision:: ident(gradVecDim,gradVecDim)
        double precision::Ut_sinvHalf_U(gradVecDim,gradVecDim)
        double precision::scfvecCheck(gradVecDim,gradVecDim)

        nbas=gradVecDim
        nbas2=int(gradVecDim**2)

        print*,'@rotateMOs_AMSgrad: rotation vec Kappa:'
        call output(kappa,1,nbas,1,nbas,nbas,nbas,1)
        ! initialize Umat with identity
        Umat=0.0d0
        do i=1,nbas
          Umat(i,i)=1.0d0
        enddo
        Umat=Umat+kappa+0.5d0*matmul(kappa,kappa)
       ! call GramSchmidt(Umat,gradVecDim,gradVecDim)

        print*,'@rotateMOs_AMSgrad: MO rotation matrix, Umat'
        call output(Umat,1,nbas,1,nbas,nbas,nbas,1)

        ! First, operate on SCFEVECA like in Brueckner
        ! Note: w/o orthogonalization like in 
        ! /vcc/zackRott1.F90
        call getrec(20,'JOBARC','SCFEVECA',NBAS2,scfVecs)
        call getrec(20,'JOBARC','SCFEVCA0',NBAS2,scfVeccheck)
        print*,'@rotateMOs_AMSgrad: Check diff b/t pre/post HF vecs...'
        call checksum('@rotateMOs: checksum HF vecs',scfVeccheck,nbas2)
        call checksum('@rotateMOs: checksum post-HF vecs',scfVecs,nbas2)

        print*,'@rotateMOs: post SCF vecs:'
        call output(scfVecs,1,nbas,1,nbas,nbas,nbas,1)

        print*,'@rotateMOs: pre SCF vecs:'
        call output(scfVeccheck,1,nbas,1,nbas,nbas,nbas,1)

        
        rotatedVecs=0.0d0
        rotatedVecs=matmul(scfVecs,Umat)
        call checksum('@rotateMO: checksumrotatedVec',rotatedVecs,nbas2)
        print*,'@rotateMOs: Printing rotated HF vecs...'
        call output(rotatedVecs,1,nbas,1,nbas,nbas,nbas,1)
        call putrec(20,'JOBARC','SCFEVECA',nbas2,rotatedVecs)

        call getrec(20,'JOBARC','AOOVRLAP',nbas2,AOovrlap)
        call AO2MO2(AOovrlap,AOovrlap,rotatedVecs,scr,nbas,nbas,1)
        call eig(AOovrlap,rotatedVecs,nbas,nbas,0)
        call invsqt(AOovrlap,nbas+1,nbas)
        Ut_sinvHalf_U=matmul(matmul(rotatedVecs,AOovrlap),&
     &                         transpose(rotatedVecs))

        ! Now transform coeffs. from MO basis to AO basis:
        call getrec(20,'JOBARC','SCFEVECA',NBAS2,finalSCFvec)
        transVec=0.0d0
        transVec=matmul(finalSCFvec,Ut_sinvHalf_U)
        call checksum('@zacRott1:chksum transVec putrec',transVec,nbas2)
        call putrec(20,'JOBARC','SCFEVCA0',nbas2,transVec)

        ! check if C^t S C ==1
        call getrec(20,'JOBARC','AOOVRLAP',nbas2,AOovrlap)
        ident=matmul(matmul(transpose(transVec),AOovrlap),transVec)
        print*,'@rotateMOs: Checking C^t S C ==1'
        call output(ident,1,nbas,1,nbas,nbas,nbas,1)
#ifdef _NO_SKIP        
        ! Load AO ovrlap matrix and transform to MO basis
        CALL GETREC(20,'JOBARC','AOOVRLAP',NBAS2,AOovrlap)
        call checksum('@rotateMO: checksum AO overlap',AOovrlap,nbas2)
        call AO2MO2(AOovrlap,AOovrlap,rotatedVecs,scr,nbas,nbas,1)
        CALL CHECKSUM('@rotateMO:Ctrans.Covrlap',AOovrlap,nbas2)
      call eig(AOovrlap,rotatedVecs,nbas,nbas,0)
      call invsqt(AOovrlap,nbas+1,nbas)
      call checksum('@rotateMO:chksum eig vecs:',rotatedVecs,nbas2)
      call checksum('@rotateMO:chksum eig vals:',AOovrlap,nbas2)
      xx=0.0d0
      xx=matmul(matmul(rotatedVecs,AOovrlap),transpose(rotatedVecs))
      call checksum('@zacRott1: chksum otho.trans.',xx,nbas2)

        ! NOW TRANSFORM MO VECTOR TO ORTHOGONAL REPRESENTATION AND
        ! DEPOSIT
        call getrec(20,'JOBARC','SCFEVECA',NBAS2,finalSCFvec)
        call checksum('@rotateMO: checksumrotatedVec',finalSCFvec,nbas2)
        transVec=0.0d0
        transVec=matmul(finalSCFvec,xx)
      call checksum('@zacRott1:chksum transVec putrec',transVec,nbas2)
      call putrec(20,'JOBARC','SCFEVCA0',nbas2,transVec)

        ! CHECK C^T S C ==1
                      
        AOovrlap=0.0d0
        CALL GETREC(20,'JOBARC','AOOVRLAP',NBAS2,AOovrlap)
        ident=matmul(matmul(transpose(transVec),AOovrlap),transVec)
        print*,'@rotateMO: check for C^tSC==1'
        call output(ident,1,nbas,1,nbas,nbas,nbas,1)
#endif
#ifdef _NO_SKIP
        print*,'@rotateMOs: Umat before GS:'
        call output(Umat,1,gradVecDim,1,gradVecDim,
     &                  gradVecDim,gradVecDim,1)
        ! Orthogonalize Umat
        call GramSchmidt(Umat,gradVecDim,gradVecDim)

        print*,'@rotateMOs: check on unitary condition Umat'
        checkUnitary=matmul(Umat,transpose(Umat))
        call output(checkUnitary,1,gradVecDim,1,gradVecDim,
     &                  gradVecDim,gradVecDim,1)        

        print*,'@rotateMOs: Umat after GS:'
        call output(Umat,1,gradVecDim,1,gradVecDim,
     &                  gradVecDim,gradVecDim,1)

        ! Load old MO vectors
        call Getrec(20,"JOBARC","SCFEVCA0",Isqlen,oldMO)


        ! call external routine to build the new MO
        ! vectors using oldMOs and Umat
        print*,'@rotateMOs: SCF vecs b4 rotation'
        call output(oldMO,1,gradVecDim,1,gradVecDim,
     &                  gradVecDim,gradVecDim,1)

        call buildNewVecs(oldMO,Umat,newMO,gradVecDim)

        print*,'@rotateMOs: SCF vecs after rotation'
        call output(newMO,1,gradVecDim,1,gradVecDim,
     &                  gradVecDim,gradVecDim,1)

        ! Finally, put the new vectors into JOBARC

        call Putrec(20,"JOBARC","SCFEVCA0",Isqlen,newMO)
#endif
        end subroutine
