      Subroutine form_h1bar(Fmi_a,Fmi_b,Fae_a,Fae_b,Fme_a,Fme_b,Fem_a,
     +                      Fem_b,H1bar_aa,H1bar_bb,Fockoo_a,Fockoo_b,
     +                      Fockvv_a,Fockvv_b,Work,Nocc_a,Nocc_b,
     +                      Nvrt_a,Nvrt_b,Nbasis,Memleft)

      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Memleft)

      Dimension Fmi_a(Nocc_a,Nocc_a)
      Dimension Fmi_b(Nocc_b,Nocc_b)
      Dimension Fae_a(Nvrt_a,Nvrt_a)
      Dimension Fae_b(Nvrt_b,Nvrt_b)
      Dimension Fme_a(Nocc_a,Nvrt_a)
      Dimension Fme_b(Nocc_b,Nvrt_b)
      Dimension Fem_a(Nvrt_a,Nocc_a)
      Dimension Fem_b(Nvrt_b,Nocc_b)

      Dimension H1bar_aa(Nocc_a*Nvrt_a,Nocc_a*Nvrt_a)
      Dimension H1bar_bb(Nocc_a*Nvrt_a,Nocc_a*Nvrt_a)

      Dimension fockoo_a(Nocc_a,Nocc_a)
      Dimension fockoo_b(Nocc_b,Nocc_b)
      Dimension fockvv_a(Nvrt_a,Nvrt_a)
      Dimension fockvv_b(Nvrt_b,Nvrt_b)

      Integer A,B,I,J

      Data Ione,Done,Dzro,Thres /1,1.0D0,0.0D0,1.0D-08/

#include "cse_vars.com"

      Ldim_a = Nocc_a*Nvrt_a
      Ldim_b = Nocc_b*Nvrt_b

      Call Dzero(H1bar_aa,Ldim_a*Ldim_a)
      Call Dzero(H1bar_bb,Ldim_b*Ldim_b)

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The Fmi_a"
      call output(Fmi_a,1,Nocc_a,1,Nocc_a,Nocc_a,Nocc_a,1)
      Write(6,"(a)") "The Fae_a"
      call output(Fae_a,1,Nvrt_a,1,Nvrt_a,Nvrt_a,Nvrt_a,1)
      Write(6,"(a)") "The Fme_a"
      call output(Fme_a,1,Nocc_a,1,Nvrt_a,Nocc_a,Nvrt_a,1)
      Write(6,"(a)") "The Fem_a"
      call output(Fme_a,1,Nvrt_a,1,Nocc_a,Nvrt_a,Nocc_a,1)
#endif

C Complete the F(mi) and F(ae) H1bar elements with the diagonal
C elements.

C OO-AA block

      Do I = 1, Nocc_a
      Do J = 1, Nocc_a
         If (I .Eq. J) Then
            E_i = Ocn_oa(I)
            Fmi_a(I,I) =  Fockoo_a(I,I)*E_i + Fmi_a(I,I)
         Else
            Fmi_a(J,I) =  Fmi_a(J,I)
         Endif
      Enddo
      Enddo

C VV-AA block

      Do A = 1, Nvrt_a
      Do B = 1, Nvrt_a
         If (A .EQ. B) Then
            E_a = (1.0D0-Ocn_va(A))
            Fae_a(A,A) = Fockvv_a(A,A)*E_a + Fae_a(A,A)
         Else
            Fae_a(B,A) = Fae_a(B,A)
         Endif
      Enddo
      Enddo

C OO-BB block

      Do I = 1, Nocc_b
      Do J = 1, Nocc_b
         If (I .Eq. J) Then
            E_i = Ocn_ob(I)
            Fmi_b(I,I) =  Fockoo_b(I,I)*E_i + Fmi_b(I,I)
         Else
            Fmi_b(J,I) =  Fmi_b(J,I)
         Endif
      Enddo
      Enddo

C VV-BB block

      Do A = 1, Nvrt_b
      Do B = 1, Nvrt_b
         If (A .EQ. B) Then
            E_a = (1.0D0-Ocn_vb(A))
            Fae_b(A,A) = Fockvv_b(A,A)*E_a + Fae_b(A,A)
         Else
            Fae_b(B,A) = Fae_b(B,A)
         Endif
      Enddo
      Enddo

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "Fmi_aa + Diagonals"
      Call output(Fmi_a,1,Nocc_a,1,Nocc_a,Nocc_a,Nocc_a,1)
      Write(6,"(a)") "Fae_aa + Diagonals"
      Call output(Fae_a,1,Nvrt_a,1,Nvrt_a,Nvrt_a,Nvrt_a,1)
      Write(6,"(a)") "Fmi_bb + Diagonals"
      Call output(Fmi_b,1,Nocc_b,1,Nocc_b,Nocc_b,Nocc_b,1)
      Write(6,"(a)") "Fae_bb + Diagonals"
      Call output(Fae_b,1,Nvrt_b,1,Nvrt_b,Nvrt_b,Nvrt_b,1)
#endif 

      Call Built_jacobi(H1bar_aa,Fmi_a,Fae_a,Nocc_a,Nvrt_a)
      Call Built_jacobi(H1bar_bb,Fmi_b,Fae_b,Nocc_b,Nvrt_b)

c Invert the H1bar matrices 

      Nrow    = Nvrt_a*Nocc_a
      Ncol    = Nvrt_a*Nocc_a

      Ieval_r = Ione
      Ieval_l = Ieval_r + Nrow
      Ievec_r = Ieval_l + Nrow 
      Ievec_l = Ievec_r + Nrow*Ncol
      Iscr1   = Ievec_l + Nrow*Ncol
      Iend    = Iscr1   + Nrow*Ncol

#ifdef _DEBUG_LVL0
      Call Dcopy(Nrow*Ncol,H1bar_aa,1,Work(Iend),1)
#endif 
   
      If (Iend .Ge. Memleft) Call Insmem("form_h1bar",Iend,Memleft)

C Alpha-alpha one electron Hamiltonian 

      Call Dgeev("V","V",Nrow,H1bar_aa,Nrow,Work(Ieval_r),Work(Ieval_l),
     +            Work(Ievec_r),Nrow,Work(Ievec_l),Nrow,Work(Iscr1),
     +            4*Nrow,Ierr) 

      Call Set_lr_norm(Work(Ievec_l),Work(Ievec_r),Work(Iscr1),
     +                 Nrow,Ncol) 

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The right eigenvalues if H1bar_aa"
      Write(6,"(6(1x,F10.5))") (Work(Ieval_r+I-1),i=1,Nrow) 
#endif 
      Call Dzero(H1bar_aa,Nrow*Ncol)
      Do I = 1, Nrow
         If (Abs(Work(Ieval_r+I-1)) .Gt. Thres) Then
             H1bar_aa(I,I) = 1.0d0/Work(Ieval_r+I-1) 
         Else
             H1bar_aa(I,I) = 0.0D0
         Endif 
      Enddo
 
      Call Dgemm("N","N",Nrow,Ncol,Ncol,Done,Work(Ievec_l),Nrow,
     +            H1bar_aa,Nrow,Dzro,Work(Iscr1),Nrow)

      Call Dgemm("N","T",Nrow,Ncol,Ncol,Done,Work(Iscr1),Nrow,
     +            Work(Ievec_r),Nrow,Dzro,H1bar_aa,Nrow)

      Nrow    = Nvrt_b*Nocc_b
      Ncol    = Nvrt_b*Nocc_b

      Ieval_r = Ione
      Ieval_l = Ieval_r + Nrow
      Ievec_r = Ieval_l + Nrow
      Ievec_l = Ievec_r + Nrow*Ncol
      Iscr1   = Ievec_l + Nrow*Ncol
      Iend    = Iscr1   + Nrow*Ncol

C#ifdef _DEBUG_LVLM
CC Invoke to heck the inverse of H1bar
C      Call Dcopy(Nrow*Ncol,H1bar_bb,1,Work(Iend),1)
C#endif
C  
C      If (Iend .Ge. Memleft) Call Insmem("form_h1bar",Iend,Memleft)

C Beta-beta one electron Hamiltonian

      Call Dgeev("V","V",Nrow,H1bar_bb,Nrow,Work(Ieval_r),Work(Ieval_l),
     +            Work(Ievec_r),Nrow,Work(Ievec_l),Nrow,Work(Iscr1),
     +            4*Nrow,Ierr)

      Call Set_lr_norm(Work(Ievec_l),Work(Ievec_r),Work(Iscr1),
     +                 Nrow,Ncol)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The right eigenvalues if H1bar_bb"
      Write(6,"(6(1x,F10.5))") (Work(Ieval_r+I-1),i=1,Nrow)
#endif

      Call Dzero(H1bar_bb,Nrow*Ncol)
      Do I = 1, Nrow
         If (Abs(Work(Ieval_r+I-1)) .Gt. Thres) Then
             H1bar_bb(I,I) = 1.0d0/Work(Ieval_r+I-1)
         Else
             H1bar_bb(I,I) = 0.0D0
         Endif
      Enddo

      Call Dgemm("N","N",Nrow,Ncol,Ncol,Done,Work(Ievec_l),Nrow,
     +            H1bar_bb,Nrow,Dzro,Work(Iscr1),Nrow)

      Call Dgemm("N","T",Nrow,Ncol,Ncol,Done,Work(Iscr1),Nrow,
     +            Work(Ievec_r),Nrow,Dzro,H1bar_bb,Nrow)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The alpha H1bar inverse matrices"
      Call output(H1bar_aa,1,Nrow,1,Ncol,Nrow,Ncol,1)
      Write(6,"(a)") "The beta  H1bar inverse matrices"
      Call output(H1bar_bb,1,Nrow,1,Ncol,Nrow,Ncol,1)

#ifdef _DEBUG_LVL0
C Invoke to check the inverse of H1bar
      Write(6,"(a)") "The (h1bar)(h1bar)^-1 "
      Call Dgemm("N","N",Nrow,Ncol,Ncol,Done,Work(Iend),Nrow,
     +            H1bar_aa,Nrow,Dzro,Work(Iscr1),Nrow)
      Call output(Work(Iscr1),1,Nrow,1,Ncol,Nrow,Ncol,1)
#endif 
#endif 

      Return
      End 
