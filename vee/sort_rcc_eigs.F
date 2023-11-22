      Subroutine Sort_rcc_eigs(Evals,Evecr,Evecl,Matdim)

      Implicit Double Precision(A-H,O-Z)
    
      Dimension Evals(Matdim),Evecr(Matdim,Matdim)
      Dimension Evecl(Matdim,Matdim)

C Sort the eigenvalues and vectors in ascending order

      Do Idim = 1, Matdim
      Do Jdim = Idim+1, Matdim

         If (Evals(Idim) .Gt. Evals(Jdim)) Then
            Eval_tmp       = Evals(Idim)
            Evals(Idim) = Evals(Jdim)
            Evals(Jdim) = Eval_tmp
            
            Do Kdim = 1, Matdim
               Evec_tmpl =  Evecl(Kdim,Idim)
               Evec_tmpr =  Evecr(Kdim,Idim)
               Evecl(Kdim,Idim) = Evecl(Kdim,Jdim)
               Evecr(Kdim,Idim) = Evecr(Kdim,Jdim)
               Evecl(Kdim,Jdim) = Evec_tmpl
               Evecr(Kdim,Jdim) = Evec_tmpr
            Enddo
         Endif 
       Enddo
       Enddo 

#ifdef _DEBUG_LVLM
       Write(6,*)
       Write(6,"(a)") " The sorted real eigenvalues"
       Write(6,"(5(1x,F15.8))") (Evals(i),i=1,MATDIM)
#endif
       Return
       End
