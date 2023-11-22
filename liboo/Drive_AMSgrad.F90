       subroutine Drive_AMSgrad(Work,Maxcor,Iuhf,OOmicroItr,OOmacroItr,&
     &                        nocc,nvirt,nvecDim)
        implicit none
        integer Maxcor
      double precision Work(Maxcor)
      integer iuhf,maxitr,OOmicroItr,OOmacroItr,nocc,nvirt,nvecDim
      double precision TOL,gradnorm,kappaOld((nocc+nvirt)**2)
      double precision kappa((nocc+nvirt)**2)
      double precision globalGradNorm,diff
      double precision gk((nocc+nvirt)**2),yk_1((nocc+nvirt)**2)
      double precision sk_1((nocc+nvirt)**2)
      double precision step,gkOld(nvecDim)
      double precision::m((nocc+nvirt)**2),v((nocc+nvirt)**2)
      double precision::vhat((nocc+nvirt)**2)
      double precision::const,eta,beta1,beta2,kappa_prev,v_prev
      double precision::gk_RMS, gk_MAX_elem,kappa_RMS,kappa_MAX_elem
      double precision::history((nocc+nvirt)**2,4),m_prev,v_hat_prev
      Parameter(eta=0.001d0)
      Parameter(beta1=0.9d0)
      Parameter(beta2=0.999d0)
      Parameter(const=10E-8)
      integer i,nbas,nbas2
      logical converged
        print*,'@Drive_AMSgrad: Initializing subroutine...'
        nbas=nocc+nvirt
        nbas2=nbas*nbas
        call getrec(20,"JOBARC","ORBRTGRD",nbas2,gk)
        print*,'orb grad,oomicroItr',OOmicroItr
        call output(gk,1,nbas,1,nbas,nbas,nbas,1)
        if (OOmicroItr.eq.1) then
          history=0.0d0
          !history(:,3)=-1.0d0*gk
          m=0.0d0
          v=0.0d0
          vhat=0.0d0
          kappa=0.0d0
        else
          call GETrec(20,"JOBARC","HISTORYV",nbas2*4,history)
          !call getrec(20,"JOBARC","MVECTOR ",nbas2,m)
          !call getrec(20,"JOBARC","VVECTOR ",nbas2,v)
          !call getrec(20,"JOBARC","VHATVECT",nbas2,vhat)
          !call getrec(20,"JOBARC","KAPPAVEC",nbas2,kappa)
        endif
        do i=1,nbas2
         m_prev=history(i,1)
         v_hat_prev=history(i,2)
         kappa_prev=history(i,3)
         v_prev=history(i,4)


         m(i)=beta1*m_prev+(1.0d0-beta1)*gk(i)
         v(i)=beta2*v_prev+(1.0d0 -beta2)*(gk(i)**2)
         print*,'v(i):',v(i)
         vhat(i)=max(v_hat_prev,v(i))
         step=(eta*m(i))/(sqrt(vhat(i))+const)

         kappa(i)=kappa_prev-(eta*m(i))/(sqrt(vhat(i))+const)
         print*,eta,m(i),vhat(i),const,kappa(i)
         history(i,1)=m(i)
         history(i,2)=vhat(i)
         history(i,3)=kappa_prev-step
         history(i,4)=v(i)
        enddo
        call putrec(20,"JOBARC","HISTORYV",nbas2*4,history)
       ! call putrec(20,"JOBARC","KAPPAVEC",nbas2,kappa)
       ! call putrec(20,"JOBARC","MVECTOR ",nbas2,m)
       ! CALL PUTREC(20,"JOBARC","VHATVECT",nbas2,vhat)
       ! CALL PUTREC(20,"JOBARC","VVECTOR ",nbas2,v) 
        call rotateMOs_AMSgrad(-gk,nbas,nbas,nocc,nvirt)!history(:,3),nbas,nbas,nocc,nvirt)


        call Get_RMS_MAX(gk,nbas2, gk_RMS, gk_MAX_elem)
        call Get_RMS_MAX(history(:,3),nbas2, kappa_RMS, kappa_MAX_elem)
        print*,'@Drive_AMSgrad: Orbital Gradient RMSE & Max_elem:',&
     &                          gk_RMS,gk_MAX_elem
        print*,'@Drive_AMSgrad: Rotation Kappa RMSE & Max_elem:',&
     &                          kappa_RMS,kappa_MAX_elem
        call Put_ConvergInfo(OOmicroItr,gk_RMS,gk_MAX_elem, &
     &                          kappa_RMS,kappa_MAX_elem)
#ifdef _NO_SKIP
       TOL=10**-5
      ! check for convergence of gradient first
       call getrec(20,"JOBARC","ORBRTGRD",nvecDim,gk)
       call checkGrad(gk,nvecDim,TOL,converged)
       gradnorm=dot_product(gk,gk)**(0.5d0)
       print*,'norm of the gradient vector: ', gradnorm
       if (gradnorm.lt.TOL) then
         call getrec(20,"JOBARC","GLOBNORM",1,globalGradNorm)
         diff=gradnorm-globalGradNorm
         if (diff.lt.TOL) then !! reach global convergence
           CALL putrec(20,"JOBARC","MACCONVG",1,.True.)
           CALL putrec(20,"JOBARC","MICCONVG",1,.True.)
           print*,'Reached global convergence!!!'
           print*,'Required:',OOmacroItr,'macro-iterations'
           return
         else !! only reach local convergence of gradient
           CALL putrec(20,"JOBARC","MICCONVG",1,.True.)
           OOmacroItr=OOmacroItr+1
           OOmicroItr=1
         endif
        endif




       if (OOmacroItr.eq.1 .and. OOmicroItr.eq.1) then
         alpha=1.0d0
         call putrec(20,"JOBARC","OLDGRAD",nvecDim,gk) 
         print*,'b4 do loop'
         kappaOld=0.0d0
!         do i=1,(nocc+nvirt)**2,nocc+nvirt+1
!          kappaOld(i)=0.0d0!1.0d0
!         kappaOld=-1.0d0*gk
         Call checksum("kappaOld", kappaOld, nvecDim)
         print*,'after do loop'
         call putrec(20,"JOBARC","OLDKAPPA",nvecDim,kappaOld)
        
       else
         call getrec(20,"JOBARC","OLDGRAD",nvecDim,gkold)
         yk_1=gk-gkold
         Call checksum("gk", gk, nvecDim)
         Call checksum("gkOld", gkOld, nvecDim)
         Call checksum("yk_1", yk_1, nvecDim)
         call putrec(20,"JOBARC","OLDGRAD",nvecDim,gk)

         call getrec(20,"JOBARC","OLDKAPPA",nvecDim,kappaOld)
         call getrec(20,"JOBARC","NEWKAPPA",nvecDim,kappaNew)
         sk_1=kappaNew-kappaOld
         Call checksum("kappaNew", kappaNew, nvecDim)
         Call checksum("kappaOld", kappaOld, nvecDim)
         call putrec(20,"JOBARC","OLDKAPPA",nvecDim,kappaNew)

         !alpha=dot_product(sk_1,sk_1)/dot_product(sk_1,yk_1)
         print*,'SD optimized alpha: ', alpha
       endif
       alpha=1.0d0
       
       kappaNew=0.0d0
!       kappaNew=KappaOld-alpha*gk
       call scale_SD_kappa(kappaNew,gk,nocc,nvirt,0.5d0)
       Call checksum("kappaNew", kappaNew, nvecDim)
       call putrec(20,"JOBARC","NEWKAPPA",nvecDim,kappaNew)
       print*,'building newcoeffs'
       call rotateMOs(kappaNew,nocc+nvirt,nocc+nvirt,nocc,nvirt)
       !call GetNewCoeffs(kappaNew,nocc,nvirt)
       print*,'done build newcoeffs'
       call getrec(20,"JOBARC","MICRONUM",1,OOmicroItr)
       OOmicroItr=OOmicroItr+1
       print*,'beginning microiteration: ', OOmicroItr
       call putrec(20,"JOBARC","MICRONUM",1,OOmicroItr)
#endif
       end subroutine
