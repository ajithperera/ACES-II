      Subroutine Set_lr_norm(Evec_r,Evec_l,Scr,Nrow,Ncol)

      Implicit Double Precision(A-H,O-Z)

      Dimension Evec_r(Nrow,Ncol),Evec_l(Nrow,Ncol)
      Dimension Scr(Nrow,Ncol)

      Data Done,Dzro/1.0D0,0.0D0/

      Do I = 1, Ncol 
      Do J = 1, Ncol
         DLnrn = Ddot(Nrow,Evec_r(1,I),1,Evec_l(1,J),1)
         Scr(I,J) = DLnrn
      Enddo 
      Enddo 

      Do I = 1, Ncol
         Scale_fac = 1.0D0/Scr(I,I) 
         Call Dscal(Nrow,Scale_fac,Evec_r(1,I),1)
      Enddo 

#ifdef _DEBUG_LVLM
      Do I = 1, Ncol 
      Do J = 1, Ncol
         DLnrn = Ddot(Nrow,Evec_r(1,I),1,Evec_l(1,J),1)
         Scr(I,J) = DLnrn
      Enddo 
      Enddo 
      Write(6,"(a)") " The <L|R> norm"
      call output(Scr,1,Nrow,1,Ncol,Nrow,Ncol,1)
#endif 

      Return
      End 
