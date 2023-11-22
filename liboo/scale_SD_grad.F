        subroutine scale_SD_kappa(kappa,grad,nocc,nvirt,max_step)
        integer nocc,nvirt
        double precision max_step
        double precision grad(nocc+nvirt,nocc+nvirt)
        double precision kappa(nocc+nvirt,nocc+nvirt)
        double precision MOfockMat(nocc+nvirt,nocc+nvirt)
        integer nbas
        logical converged

        nbas=nocc+nvirt
        kappa=grad
c        print*,"@scale_SD_grad: Current kappa:"
c        call output(kappa,1,nbas,1,nbas,nbas,nbas,1)
c
c        call getrec(20,"JOBARC","MOFOCKMX",(nocc+nvirt)**2,MOfockMat)
c        print*,"@scale_SD_grad: MO Fock matrix:"
c        call output(MOfockMat,1,nbas,1,nbas,nbas,nbas,1)        


        call ScaleKappaOV(kappa(1:nocc,nocc+1:nbas),nocc,nvirt,max_step,
     &                          MOfockMat)

        print*,'@scale_SD_grad: Gradient matrix before scale:'
        call output(grad,1,nbas,1,nbas,nbas,nbas,1)
        kappa(nocc+1:nbas,1:nocc)=-1.0d0*transpose(kappa(1:nocc,
     &                                                  nocc+1:nbas))
        print*,'@scale_SD_grad: Rescaled kappa'
        call output(kappa,1,nbas,1,nbas,nbas,nbas,1)
        
        ! check RMS of kappa
        print*,'@scale_SD_grad: Checking RMS of kappa:'
        call checkGrad(kappa,nbas**2,10**-5,converged)

        end subroutine

