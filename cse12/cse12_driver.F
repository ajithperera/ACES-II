      Subroutine Cse12_driver(W0_aa,W0_bb,W0_ab,W1_aa,W1_bb,
     +                        W1_ab,W2_aa,W2_bb,W2_ab,
     +                        W3_aaaa,W3_bbbb,W3_abab,W3_baba,
     +                        W3_abba,W3_baab,W_aa,W_bb,W_ab,
     +                        W4_aa,W4_bb,W4_ab,W4_ba,
     +                        W5_aa,W5_bb,W5_ab,W5_ba,
     +                        Wmnij_1,Wmnij_2,Wmnij_3,Wabef_1,
     +                        Wabef_2,Wabef_3,Wmbej_1,Wmbej_2,
     +                        Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,
     +                        T2old_aa,T2old_bb,T2old_ab,
     +                        T1taup_aa,T1taup_bb,T1taup_ab,T1tau_aa,
     +                        T1tau_bb,T1tau_ab,T1old_aa,T1old_bb,
     +                        Fockoo_a,Fockoo_b,Fockvv_a,
     +                        Fockvv_b,Fockov_a,Fockov_b,
     +                        Fockoo_od_a,Fockoo_od_b,
     +                        Fockvv_od_a,Fockvv_od_b,Fmi_a,
     +                        Fmi_b,Fae_a,Fae_b,Fme_a,Fme_b,
     +                        Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +                        P,Work,Maxcor,Nbasis,Length,
     +                        Ref_energy,Frac_occ,Open_shell)

      Implicit Double Precision(A-H,O-Z)

#include "cse12_vars.com"

      Dimension Work(Maxcor)
      Dimension P(Length) 

      Dimension T2old_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2old_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2old_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T1taup_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T1taup_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T1taup_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T1tau_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T1tau_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T1tau_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T1old_aa(Nvrt_a,Nocc_a)
      Dimension T1old_bb(Nvrt_b,Nocc_b)

      Dimension W0_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Wc_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension W0_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W1_aa(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension W1_bb(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension W1_ab(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension W2_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension W2_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension W2_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension W3_aaaa(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W3_bbbb(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W3_abab(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W3_baba(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension W3_abba(Nocc_a,Nvrt_b,Nvrt_b,Nocc_a)
      Dimension W3_baab(Nocc_b,Nvrt_a,Nvrt_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_b)

      Dimension W5_aa(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension W5_bb(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension W5_ab(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension W5_ba(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)

      Dimension Wmnij_1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Wmnij_2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Wmnij_3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension Wabef_1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Wabef_2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Wabef_3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension Wmbej_1(Nocc_a,Nocc_a,Nvrt_a,Nocc_b)
      Dimension Wmbej_2(Nocc_b,Nocc_b,Nvrt_b,Nocc_b)
      Dimension Wmbej_3(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)
      Dimension Wmbej_4(Nocc_b,Nocc_a,Nvrt_b,Nocc_a)
      Dimension Wmbej_5(Nocc_b,Nocc_a,Nvrt_a,Nocc_b)
      Dimension Wmbej_6(Nocc_a,Nocc_b,Nvrt_b,Nocc_a)

      Dimension Fockoo_a(Nocc_a,Nocc_a)
      Dimension Fockoo_b(Nocc_b,Nocc_b)
      Dimension Fockvv_a(Nvrt_a,Nvrt_a)
      Dimension Fockvv_b(Nvrt_b,Nvrt_b)
      Dimension Fockov_a(Nocc_a,Nvrt_a)
      Dimension Fockov_b(Nocc_b,Nvrt_b)

      Dimension Fockoo_od_a(Nocc_a,Nocc_a)
      Dimension Fockoo_od_b(Nocc_b,Nocc_b)
      Dimension Fockvv_od_a(Nvrt_a,Nvrt_a)
      Dimension Fockvv_od_b(Nvrt_b,Nvrt_b)

      Dimension Fmi_a(Nocc_a,Nocc_a)
      Dimension Fmi_b(Nocc_b,Nocc_b)
      Dimension Fae_a(Nvrt_a,Nvrt_a)
      Dimension Fae_b(Nvrt_b,Nvrt_b)
      Dimension Fme_a(Nocc_a,Nvrt_a)
      Dimension Fme_b(Nocc_b,Nvrt_b)

      Logical Converged,Fme_on,T_guess_exist 
      Logical Frac_occ,Open_shell 
      Character*80 Fullname

      Call Guesst1(T1old_aa,T1old_bb,Nocc_a,Nocc_b,Nvrt_a,
     +             Nvrt_b)

      Call Guesst2(T2old_aa,T2old_bb,T2old_ab,W0_aa,W0_bb,W0_ab,
     +             Fockoo_a,Fockoo_b,Fockvv_a,Fockvv_b,
     +             Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call Energy_frmt2(T2old_aa,T2old_bb,T2old_ab,W0_aa,W0_bb,W0_ab,
     +                  Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,E2_aa,E2_bb,E2_ab)
 
      Embpt2   = (E2_aa + E2_bb)*0.25D0 + E2_ab

      Write(6,"(a,F15.10)") "The MBPT(2) energy   : ", Embpt2

C Check whether the TGUESS file exists. Then read it and use those
C as the starting guess

      Call Gfname("TGUESS",Fullname,Namelen)
      Inquire(File=Fullname(1:Namelen),Exist=T_guess_exist)

      Iunit = 94
      If (T_guess_exist .and. .Not. Brueck) Then
         Open(Unit=Iunit,File=Fullname(1:Namelen),Status="Old",
     +        Form="Unformatted",Access="Sequential")
         Call Read_amps(Iunit,T1old_aa,T1old_bb,T2old_aa,T2old_bb,
     +                  T2old_ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,Work
     +                  Maxcor)
         Close(Unit=Iunit,Status="Keep") 
      Endif

      L_aa   = Nocc_a*Nvrt_a
      L_bb   = Nocc_b*Nvrt_b
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
     
      Ihtem_aa    = 1
      Ihtem_bb    = Ihtem_aa    + L_aa
      Ihtabij_aa  = Ihtem_bb    + L_bb
      Ihtabij_bb  = Ihtabij_aa  + L_aaaa
      Ihtabij_ab  = Ihtabij_bb  + L_bbbb 
      It1resid_aa = Ihtabij_ab  + L_abab 
      It1resid_bb = It1resid_aa + L_aa 
      It2resid_aa = It1resid_bb + L_bb
      It2resid_bb = It2resid_aa + L_aaaa
      It2resid_ab = It2resid_bb + L_bbbb
      Iend        = It2resid_ab + L_abab

      
      If (Iend .GT. Maxcor) Call Insmem("ccsd_light_driver",Iend,
     +                                   Maxcor)

      Memleft = Maxcor - Iend 
      Iter      = 0
      Converged = .False. 
      E_old     = Embpt2
      E_corr(0) = Embpt2 

#ifdef _DEBUG_LVLM
      Call Dzero(T2old_aa,L_aaaa)
      Call Dzero(T2old_bb,L_bbbb)
      Call Dzero(T2old_ab,L_abab)
      E_corr(0) = 0.0D0
      E_old     = 0.0D0
#endif
      Do while ((.NOT. Converged) .AND. Iter .Le. Cc_maxcyc) 
         Iter = Iter + 1

#ifdef _DEBUG_LVLM
      Call Guesst1(T1old_aa,T1old_bb,Nocc_a,Nocc_b,Nvrt_a,
     +             Nvrt_b)
#endif 
        If (Active_space) Then
           Call Active_sp_proj(T1old_aa,T1old_bb,T2old_aa,T2old_bb,
     +                         T2old_ab,P,Length,Nocc_a,Nocc_b,
     +                         Nvrt_a,Nvrt_b) 
        Endif 

C Form Htilde (Ht=({H}e^T1)c) (One body similarity transformed 
C Hamiltonian (for legacy reason Fmi,Fae,Fme intermediate names 
C are kept in the calling routine; I may have renamed them in
C the routines called by the driver routines).

         Call Form_t1taup(T1taup_aa,T1taup_bb,T1taup_ab,T1old_aa,
     +                    T1old_bb,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

         Call Form_h1tilde(T1old_aa,T1old_bb,T1taup_aa,T1taup_bb,
     +                     T1taup_ab,Fmi_a,Fmi_b,Fae_a,Fae_b,Fme_a,
     +                     Fme_b,Work(Ihtem_aa),Work(Ihtem_bb),W_aa,
     +                     W_bb,W_ab,W3_aaaa,W3_bbbb,W3_abab,W3_baba,
     +                     W3_abba,W3_baab,W4_aa,W4_bb,W4_ab,
     +                     W4_ba,W5_aa,W5_bb,W5_ab,W5_ba,Fockoo_a,
     +                     Fockoo_b,Fockvv_a,Fockvv_b,Fockoo_od_a,
     +                     Fockoo_od_b,Fockvv_od_a,Fockvv_od_b,
     +                     Fockov_a,Fockov_b,Nocc_a,Nocc_b,Nvrt_a,
     +                     Nvrt_b)

C Form right hand side of the singles residual, R(a,i)

         Call Form_rhs_h1trai(Work(Ihtem_aa),Work(Ihtem_bb),Fme_a,Fme_b,
     +                        T1old_aa,T1old_bb,T2old_aa,T2old_bb,
     +                        T2old_ab,W_aa,W_bb,W_ab,W4_aa,W4_bb,
     +                        W4_ab,W4_ba,W5_aa,W5_bb,W5_ab,W5_ba,
     +                        Work(It1resid_aa),Work(It1resid_bb),
     +                        Fockoo_a,Fockoo_b,Fockvv_a,Fockvv_b,
     +                        Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

         Call Form_t1tau(T1tau_aa,T1tau_bb,T1tau_ab,T1old_aa,T1old_bb,
     +                   Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

         Call Form_h2tilde(T1old_aa,T1old_bb,T1tau_aa,T1tau_bb,T1tau_ab,
     +                     W0_aa,W0_bb,W0_ab,W1_aa,W1_bb,W1_ab,W_aa,
     +                     W_bb,W_ab,W2_aa,W2_bb,W2_ab,W3_aaaa,W3_bbbb,
     +                     W3_abab,W3_baba,W3_abba,W3_baab,W4_aa,W4_bb,
     +                     W4_ab,W4_ba,W5_aa,W5_bb,W5_ab,W5_ba,Wmnij_1,
     +                     Wmnij_2,Wmnij_3,Wabef_1,Wabef_2,Wabef_3,
     +                     Wmbej_1,Wmbej_2,Wmbej_3,Wmbej_4,Wmbej_5,
     +                     Wmbej_6,Work(Ihtabij_aa),Work(Ihtabij_bb),
     +                     Work(Ihtabij_ab),Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Form modified Htae and Htmi intermediates (increment with Htme)

          Call form_modf_f(Fae_a,Fae_b,Fmi_a,Fmi_b,Fme_a,Fme_b,T1old_aa,
     +                     T1old_bb,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,1.0D0)

C Form right hand side of the doubles residual from Htilde Hamiltonian
C elements, R(ij,ab)

          Call Form_rhs_h2trabij(T2old_aa,T2old_bb,T2old_ab,Fmi_a,Fmi_b,
     +                           Fae_a,Fae_b,Wmnij_1,Wmnij_2,Wmnij_3,
     +                           Wabef_1,Wabef_2,Wabef_3,Wmbej_1,
     +                           Wmbej_2,Wmbej_3,Wmbej_4,Wmbej_5,
     +                           Wmbej_6,Work(Ihtabij_aa),
     +                           Work(Ihtabij_bb),Work(Ihtabij_ab),
     +                           Work(It2resid_aa),Work(It2resid_bb),
     +                           Work(It2resid_ab),Nocc_a,Nocc_b,
     +                           Nvrt_a,Nvrt_b)
C
C Form the Fmi,Fae,Wmnij,Wabef,Wmbej again, but now the contributions
C only come from T2. I can overwrite the existing locations. 

          Call Form_t2intms(T2old_aa,T2old_bb,T2old_ab,W_aa,W_bb,W_ab,
     +                      Fmi_a,Fmi_b,Fae_a,Fae_b,Wmnij_1,
     +                      Wmnij_2,Wmnij_3,Wabef_1,Wabef_2,Wabef_3,
     +                      Wmbej_1,Wmbej_2,Wmbej_3,Wmbej_4,Wmbej_5,
     +                      Wmbej_6,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

          Call Form_rhs_t2rabij(T2old_aa,T2old_bb,T2old_ab,Fmi_a,
     +                          Fmi_b,Fae_a,Fae_b,Wmnij_1,Wmnij_2,
     +                          Wmnij_3,Wabef_1,Wabef_2,Wabef_3,
     +                          Wmbej_1,Wmbej_2,Wmbej_3,Wmbej_4,
     +                          Wmbej_5,Wmbej_6,Work(It2resid_aa),
     +                          Work(It2resid_bb),Work(It2resid_ab),
     +                          Fockoo_a,Fockoo_b,Fockvv_a,Fockvv_b,
     +                          Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

#ifdef _DEBUG
C This was how these terms were treated before these contributions
C were formulated to work T2*H_tilde. The latter form is the correct
C way to do them in Htilde formulation of CC.

         Call Form_rhs_t1t1rabij(T1tau_aa,T1tau_bb,T1tau_ab,
     +                           Wmnij_1,Wmnij_2,Wmnij_3,Wabef_1,
     +                           Wabef_2,Wabef_3,Work(It2resid_aa),
     +                           Work(It2resid_bb),Work(It2resid_ab),
     +                           Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

         Call Form_t1rai(T1old_aa,T1old_bb,Fmi_a,Fmi_b,Fae_a,Fae_b,
     +                   Work(It1resid_aa),Work(It1resid_bb),
     +                   Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
#endif 

          Call Form_residuals(T1old_aa,T1old_bb,T2old_aa,T2old_bb,
     +                        T2old_ab,Work(It1resid_aa),
     +                        Work(It1resid_bb),Work(It2resid_aa),
     +                        Work(It2resid_bb),Work(It2resid_ab),
     +                        Fockoo_a,Fockoo_b,Fockvv_a,Fockvv_b,
     +                        Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,Rms)

C T2old->T2new 

          Call Newt2(T2old_aa,T2old_bb,T2old_ab,Work(It2resid_aa),
     +               Work(It2resid_bb),Work(It2resid_ab),Fockoo_a,
     +               Fockoo_b,Fockvv_a,Fockvv_b,Nocc_a,Nocc_b,
     +               Nvrt_a,Nvrt_b)

C T1old->T1new 
          Call Newt1(T1old_aa,T1old_bb,Work(It1resid_aa),
     +               Work(It1resid_bb),Fockoo_a,Fockoo_b,Fockvv_a,
     +               Fockvv_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

          Call Energy_frmt2(T2old_aa,T2old_bb,T2old_ab,W0_aa,W0_bb,
     +                      W0_ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +                      E2ccsd_aa,E2ccsd_bb,E2ccsd_ab)

          Call Energy_frmt1(T1old_aa,T1old_bb,W0_aa,W0_bb,W0_ab,
     +                      Fockov_a,Fockov_b,Nocc_a,Nocc_b,Nvrt_a,
     +                      Nvrt_b,E1ccsd_aa,E1ccsd_bb,E1ccsd_ab,
     +                      Enhf_aa,Enhf_bb)

          Call Form_cumldens(T1old_aa,T1old_bb,Fockoo_a,Fockoo_b,
     +                       Work(Iend),Memleft,Nbasis,Ref_energy,
     +                       Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,Frac_occ)

          E2_new = (E2ccsd_aa + E2ccsd_bb)*0.25D0 + E2ccsd_ab
          E1_new = (E1ccsd_aa + E1ccsd_bb)*0.25D0 + E1ccsd_ab
          E3_new = Enhf_aa + Enhf_bb
          E_new  = E2_new + E1_new + E3_new 
          E_diff = Dabs(E_new - E_old)
          E_corr(Iter) = E_new 

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
CSSS          If (Iter .Eq. 3) Converged = .True.
          If (E_diff .LT. CC_conv) Converged = .True.
          E_old = E_new 

      enddo

      Call Getrec(20,"JOBARC","SCFENEG ",1,E_scf)
      If (Frac_occ) Then
         Write(6,"(a)") "The FON reference energy is added." 
         E_scf = Ref_energy
      Endif 

      If (Converged) Then  
         Call dump_amps(T1old_aa,T1old_bb,T2old_aa,T2old_bb,T2old_ab,
     +                  Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
         Write(6,*)
         Write(6,9) 
         Write(6,10)
         Write(6,11)
         Write(6,10)
      Else
         Write(6,*)
         Write(6,"(a,1x,I3,a)") " The CCSD did not converged in", 
     +                          cc_maxcyc, " iterations."
         Write(6,10)
      Endif 

 9    Format(T3, '    Summary of iterative solution of CCSD equations')
 10   Format(6x,59("-")) 
 11   Format(T24,'Correlation',T55,'Total',/,T8,'Iteration',T26,
     +       'Energy',T55,'Energy')
      Do I = 0, Iter
         Write(6,"(5x,I3,10x,F18.12,10x,F19.12)") I, E_corr(i), 
     +                                             E_corr(i)+E_scf
      Enddo
      Write(6,10)

      Return
      End 

