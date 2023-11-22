        subroutine checkGrad(gradVec,nvecDim,TOL,converged)
        integer nvecDim
        double precision gradVec(nvecDim),TOL,norm
        logical converged
        integer i
        norm=0.0d0!dot_product(gradVec,gradVec)**(0.5d0)
        do i=1,nvecDim
          norm=norm+gradVec(i)**2.0d0
        enddo
        norm=norm/nvecDim
        norm=(norm)**(0.5d0)
        write(6,"(a,21x,e24.12)"),'@checkGrad: RMS norm:',norm
        if (norm.lt.TOL) then
          write(6,"(a)"),'Gradient vector beneath threshold and', 
     +                   ' converged!'
          converged=.True.
        else
          converged=.False.
        endif
        end subroutine
