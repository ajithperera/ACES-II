        subroutine ScaleKappaOV(kappa,grad,nocc,nvirt,max_step)
        implicit none
        integer nocc,nvirt
        double precision max_step
        double precision kappa(nvirt,nocc)
        double precision grad(nvirt,nocc)
        double precision evals(nocc+nvirt)
        double precision denomShift,max_kappa
        integer i,occidx,vrtidx,nbas,a
        ! **NOTE** THIS SUBROUTINE ONLY RESCALES THE VO PORTION OF THE
        ! GRADIENT IN ORDER TO FORM THE VO PART OF ROTATION MATRIX
        ! Inspired by DOI: 10.1063/1.4816628

        nbas=nocc+nvirt
        ! Harvest Fock matrix eigenvalues to rescale gradient
        call getrec(20,"JOBARC","SCFEVLA0",nbas,evals)
        call checksum('@ScaleKappaOV:scf evals:    ',evals,nbas)
        call checksum('@ScaleKappaOV: grad    :    ',grad,nbas)

        ! Rescale: kappa=-grad/(2*[faa - fii])
        do i=1,nocc
          do a=1,nvirt
           denomShift=-evals(i)+evals(nocc+a)
CSSS           print*,'Denom shift: ', denomShift
           kappa(a,i)=-1.0d0*grad(a,i)
CSSS/(2.0d0*denomShift)
        enddo
        enddo
        call checksum('@ScaleKappaOV: kappa_ai:     ',kappa,nbas)

        ! Logic to rescale according to max Kappa & max_step
        ! Copied from Psi4 source, but commented out because it is not
        ! relevant for our current small problems.


!        max_kappa=0.0d0
!        do i=1,nocc*nvirt
!          if (abs(kappa(i)).gt.max_kappa) then
!                max_kappa=kappa(i)
!          endif
!        enddo
!        print*,'@ScaleKappaOV: maximum element of kappa:',max_kappa     
!        if (abs(max_kappa).gt.max_step) then
!          do i=1,nocc*nvirt
!            kappa(i)=kappa(i)*(max_step/max_kappa)
!          enddo
!        endif
        end subroutine
