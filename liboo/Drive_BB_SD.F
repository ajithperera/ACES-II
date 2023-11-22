       subroutine Drive_BB_SD(Work,Maxcor,Iuhf,OOmicroItr,OOmacroItr,
     +                        nocc,nvirt,nvecDim)
      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxcor)
      integer maxitr,OOmicroItr,OOmacroItr,nocc,nvirt,nvecDim
      double precision TOL,gradnorm,kappaOld((nocc+nvirt)**2)
      double precision kappaNew((nocc+nvirt)**2)
      double precision globalGradNorm,diff
      double precision gk((nocc+nvirt)**2),yk_1((nocc+nvirt)**2)
      double precision sk_1((nocc+nvirt)**2)
      double precision gkOld(nvecDim)
      integer i
      logical converged

       TOL=10**-5
      ! check for convergence of gradient first
       call getrec(20,"JOBARC","ORBRTGRD",nvecDim,gk)
       call checkGrad(gk,nvecDim,TOL,converged)
       gradnorm=dot_product(gk,gk)**(0.5d0)
       write(6,"(2a,e24.12)") '@-Drve_BB_SD norm of the gradient',
     +                        ' vector: ',gradnorm
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

        ! call routine that loads gradient, builds kappa,
        ! rotates/stores vectors
        call zackTESTrott1(nocc,nvirt)

       call getrec(20,"JOBARC","MICRONUM",1,OOmicroItr)
       OOmicroItr=OOmicroItr+1
       write(6,"(a,i3)"),'@-Drve_BB_SD beginning microiteration: ',
     +                     OOmicroItr
       call putrec(20,"JOBARC","MICRONUM",1,OOmicroItr)


       end subroutine
