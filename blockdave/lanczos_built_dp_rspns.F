      Subroutine Lanczos_built_dp_respns(Eig_vals,Leig_vecs,Reig_vecs,
     +                                Residual,Tm_k,Scr,Nrow,Ncol,
     +                                Irrepx,Nsize,Iuhf,Memleft)

      Implicit Double Precision (A-H,O-Z)
      Double Precision Leig_vecs
      Double Precision Null

      Dimension Eig_vals(Nrow),Leig_vecs(Nrow,Ncol),Scr(Memleft),
     +          Reig_vecs(Nrow,Ncol),Tm_k(Nrow),Residual(Nrow)

      Common /Init_norms/Sqrt_Norm,Rl_norm
      Data Autoev /27.21138602/

#ifdef _DEBUG_LVL0
      Write(6,"(a,2F15.8)") "@-Entry; Sqrt_norm,RL_Norm: ", 
     +                       Sqrt_Norm,Rl_norm
#endif
 
      Do K = 1, Ncol
         Tm_k(k) = Leig_vecs(K,1) * Reig_vecs(1,K)
         Tm_k(K) = Sqrt_Norm * Rl_norm * Tm_k(K)
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,"(a)") "The transition moments"
      Write(6,"(5(1x,F15.7))") (Tm_k(I),I=1,Ncol)
#endif
C Reorder the eigenvalues ascending order along with the transition moments
C and residuals.

      Do I = 1, Ncol
         Do J = I+1, Ncol
            If (Eig_vals(J) .LE. Eig_vals(I)) Then
C Eigenvalues 
                Eig_tmp = Eig_vals(I)
                EIg_vals(I) = Eig_vals(J) 
                Eig_vals(J) = Eig_tmp

C Transition moments 
                Tran_m_tmp = Tm_k(I)
                Tm_k(I)    = Tm_k(J)
                Tm_k(J)    = Tran_m_tmp
C Residual 
                Res_tmp    = Residual(I)
                Residual(I)= Residual(J)
                Residual(J)= Res_tmp
            Endif 
         Enddo
      Enddo

C Make a tables of these reordered values.

      Write(6,*)
      Write(6,"(4a)") "   Eigenvalues (a.u.) "," Eigenvalues (eV)",
     +                 " Transition moments","   Residual Norms"
      Write(6,"(a,a)") "   --------------------------------------",
     +                 "---------------------------------"
      Do I = 1, Ncol
         Write(6, "(2(1x,F15.7),2(5x,F15.7))") Eig_vals(I),
     +                                      Eig_vals(I)*Autoev,
     +                                      Tm_k(I),
     +                                      Residual(I)
      Enddo
      Write(6,"(a,a)") "   --------------------------------------",
     +                 "---------------------------------"

      Return
      End 
     
