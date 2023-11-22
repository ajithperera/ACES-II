      Subroutine Mbpt2_brueckner(T2aa,T2bb,T2ab,T1aa,T1bb,W_aa,W_bb,
     +                           W_ab,W0_aa,W0_bb,W0_ab,W4_aa,W4_bb,
     +                           W4_ab,W4_ba,W5_aa,W5_bb,W5_ab,W5_ba,
     +                           Fockoo_a,Fockoo_b,Fockvv_a,Fockvv_b,
     +                           Fockov_a,Fockov_b,Fme_a,Fme_b,Embpt2,
     +                           Work,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +                           Maxcor)

      Implicit Double Precision(A-H,O-Z) 

#include "ccsdlight_vars.com"

      Dimension Work(Maxcor)
      Logical T1_converged,Fme_on

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension W0_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension w0_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension W0_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_b)

      Dimension W5_aa(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension W5_bb(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension W5_ab(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension W5_ba(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension fockoo_a(Nocc_a,Nocc_a)
      Dimension fockoo_b(Nocc_b,Nocc_b)
      Dimension fockvv_a(Nvrt_a,Nvrt_a)
      Dimension fockvv_b(Nvrt_b,Nvrt_b)
      Dimension fockov_a(Nocc_a,Nvrt_a)
      Dimension fockov_b(Nocc_b,Nvrt_b)

      Dimension fme_a(Nocc_a,Nvrt_a)
      Dimension fme_b(Nocc_b,Nvrt_b)

      L_aa   = Nocc_a*Nvrt_a
      L_bb   = Nocc_b*Nvrt_b
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b

      It1resid_aa = 1
      It1resid_bb = It1resid_aa + L_aa
      It2resid_aa = It1resid_bb + L_bb
      It2resid_bb = It2resid_aa + L_aaaa
      It2resid_ab = It2resid_bb + L_bbbb
      Iend        = It2resid_ab + L_abab

      If (Iend .GT. Maxcor) Call Insmem("mbpt2_brueckner",Iend,
     +                                   Maxcor)

      Iter         = 0
      E_old        = Embpt2
      E_corr(0)    = Embpt2
      T1_converged = .False.
      
      Do while ((.NOT. T1_converged) .AND. Iter .Le. Cc_maxcyc)
         Iter = Iter + 1

C Fme(m,e) intermediate

         Call Form_fme(Fme_a,Fme_b,T1aa,T1bb,W_aa,W_bb,
     +                 W_ab,Fockov_a,Fockov_b,Nocc_a,Nocc_b,
     +                 Nvrt_a,Nvrt_b,Fme_on)
        
C Initialize the to <ab||ij>  and f(a,i)

         Call init_2abij(Work(It2resid_aa),Work(It2resid_bb),
     +                   Work(It2resid_ab),Work(It1resid_aa),
     +                   Work(It1resid_bb),W0_aa,W0_bb,W0_ab,
     +                   Fockov_a,Fockov_b,Nocc_a,Nocc_b,
     +                   Nvrt_a,Nvrt_b)

C Fme->T1

         Call Fme_int1(T2aa,T2bb,T2ab,Work(It1resid_aa),
     +                 Work(It1resid_bb),Fme_a,Fme_b,Nocc_a,Nocc_b,
     +                 Nvrt_a,Nvrt_b)

C T2->T1

         Call T2_int1(T2aa,T2bb,T2ab,W4_aa,W4_bb,W4_ab,W4_ba,W5_aa,
     +                W5_bb,W5_ab,W5_ba,Work(It1resid_aa),
     +                Work(It1resid_bb),Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C T2old->T2new

         Call Newt2(T2aa,T2bb,T2ab,Work(It2resid_aa),
     +              Work(It2resid_bb),Work(It2resid_ab),Fockoo_a,
     +              Fockoo_b,Fockvv_a,Fockvv_b,Nocc_a,Nocc_b,
     +              Nvrt_a,Nvrt_b)

C T1old->T2new

         Call Newt1(T1aa,T1bb,Work(It1resid_aa),
     +              Work(It1resid_bb),Fockoo_a,Fockoo_b,Fockvv_a,
     +              Fockvv_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
 
         Call Energy_frmt2(T2aa,T2bb,T2ab,W0_aa,W0_bb,
     +                     W0_ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +                     E2ccsd_aa,E2ccsd_bb,E2ccsd_ab)
 
         Call Energy_frmt1(T1aa,T1bb,W0_aa,W0_bb,W0_ab,Fockov_a,
     +                     Fockov_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +                     E1ccsd_aa,E1ccsd_bb,E1ccsd_ab)

         E2_new = (E2ccsd_aa + E2ccsd_bb)*0.25D0 + E2ccsd_ab
         E1_new = (E1ccsd_aa + E1ccsd_bb)*0.25D0 + E1ccsd_ab
         E_new  = E2_new + E1_new

#ifdef _DEBUG_LVL0
         Write(6,"(a,3(1x,F15.10))")"E1p2aa,E1p2bb,E1p2ab :",
     +                               (E2ccsd_aa+E1ccsd_aa)*0.25D0,
     +                               (E2ccsd_bb+E1ccsd_bb)*0.25D0,
     +                                E2ccsd_ab+E1ccsd_ab
         Write(6,"(a,2(1x,F15.10),2x,a)")
     +         "The energy and the energy difference: ", E_new, E_diff,
     +          "at"
         Write(6,"(a,1x,I3))") "iteration no.: ", Iter
#endif
         E_diff= Dabs(E_new - E_old)
         E_corr(Iter) = E_new
         IF (E_diff .LT. CC_conv) T1_converged = .True.

         E_old = E_new

      Enddo 

      If (T1_converged) Then
         Write(6,*)
         Write(6,9)
         Write(6,10)
         Write(6,11)
         Write(6,10)
      Else
         Write(6,*)
         Write(6,10)
         Write(6,"(a,1x,I3,a)") " The B-MBPT(2) did not converged in",
     +                          cc_maxcyc, " iterations."
      Endif

 9    Format(8X,"Summary of iterative solution of B-MBPT(2) equations")
 10   Format(6x,59("-"))
 11   Format(T24,'Correlation',T55,'Total',/,T8,'Iteration',T26,
     +       'Energy',T55,'Energy')
      Do I = 0, Iter
         Write(6,"(5x,I3,10x,F18.12,10x,F19.12)") I, E_corr(i),
     +                                             E_corr(i)+E_scf
      Enddo
      Write(6,10)

      Call Mbpt2_brueckner_rot(T1aa,T1bb,Work,Nocc_a,Nocc_b,Nvrt_a,
     +                         Nvrt_b,Maxcor)

      Return
      End 

