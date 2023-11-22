      Subroutine Ccsd_light_driver(W0_aa,W0_bb,W0_ab,W1_aa,W1_bb,
     +                             W1_ab,W2_aa,W2_bb,W2_ab,
     +                             W3_aaaa,W3_bbbb,W3_abab,W3_baba,
     +                             W3_abba,W3_baab,W_aa,W_bb,W_ab,
     +                             W4_aa,W4_bb,W4_ab,W4_ba,
     +                             W5_aa,W5_bb,W5_ab,W5_ba,
     +                             Wmnij_1,Wmnij_2,Wmnij_3,Wabef_1,
     +                             Wabef_2,Wabef_3,Wmbej_1,Wmbej_2,
     +                             Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,
     +                             Wabci_1,Wabci_2,Wabci_3,Wabci_4,
     +                             Waibc_1,Waibc_2,Waibc_3,Waibc_4,
     +                             Wijka_1,Wijka_2,Wijka_3,Wijka_4,
     +                             Wiajk_1,Wiajk_2,Wiajk_3,Wiajk_4,
     +                             Wijab_1,Wijab_2,Wijab_3,
     +                             T2old_aa,T2old_bb,T2old_ab,
     +                             Tau_aa,Tau_bb,Tau_ab,Tap_aa,
     +                             Tap_bb,Tap_ab,T1old_aa,T1old_bb,
     +                             Fockoo_a,Fockoo_b,Fockvv_a,
     +                             Fockvv_b,Fockov_a,Fockov_b,
     +                             Fockoo_od_a,Fockoo_od_b,
     +                             Fockvv_od_a,Fockvv_od_b,Fmi_a,
     +                             Fmi_b,Fae_a,Fae_b,Fme_a,Fme_b,
     +                             Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +                             P,Work,Maxcor,Nbasis,Length,
     +                             Ref_energy,Frac_occ,Open_shell)

      Implicit Double Precision(A-H,O-Z)

#include "ccsdlight_vars.com"

      Dimension Work(Maxcor)
      Dimension P(Length) 

      Dimension T2old_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2old_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2old_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Tau_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Tau_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Tau_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Tap_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Tap_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Tap_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

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

      Dimension Wijab_1(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension Wijab_2(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension Wijab_3(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension Wijka_1(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension Wijka_2(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension Wijka_3(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension Wijka_4(Nocc_b,Nocc_a,Nocc_b,Nvrt_a)

      Dimension Wiajk_1(Nocc_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Wiajk_2(Nocc_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Wiajk_3(Nocc_a,Nvrt_b,Nocc_a,Nocc_b)
      Dimension Wiajk_4(Nocc_b,Nvrt_a,Nocc_b,Nocc_a)

      Dimension Wabci_1(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Wabci_2(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Wabci_3(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Wabci_4(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
   
      Dimension Waibc_1(Nvrt_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension Waibc_2(Nvrt_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension Waibc_3(Nvrt_a,Nocc_b,Nvrt_a,Nvrt_b)
      Dimension Waibc_4(Nvrt_b,Nocc_a,Nvrt_b,Nvrt_a)

      Dimension Wmbej_1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Wmbej_2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Wmbej_3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Wmbej_4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension Wmbej_5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension Wmbej_6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Dimension Wijkl_1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Wijkl_2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Wijkl_3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension fockoo_a(Nocc_a,Nocc_a)
      Dimension fockoo_b(Nocc_b,Nocc_b)
      Dimension fockvv_a(Nvrt_a,Nvrt_a)
      Dimension fockvv_b(Nvrt_b,Nvrt_b)
      Dimension fockov_a(Nocc_a,Nvrt_a)
      Dimension fockov_b(Nocc_b,Nvrt_b)

      Dimension fockoo_od_a(Nocc_a,Nocc_a)
      Dimension fockoo_od_b(Nocc_b,Nocc_b)
      Dimension fockvv_od_a(Nvrt_a,Nvrt_a)
      Dimension fockvv_od_b(Nvrt_b,Nvrt_b)

      Dimension fmi_a(Nocc_a,Nocc_a)
      Dimension fmi_b(Nocc_b,Nocc_b)
      Dimension fae_a(Nvrt_a,Nvrt_a)
      Dimension fae_b(Nvrt_b,Nvrt_b)
      Dimension fme_a(Nocc_a,Nvrt_a)
      Dimension fme_b(Nocc_b,Nvrt_b)

      Logical Converged,Fme_on,T_guess_exist 
      Logical Frac_occ,Open_shell
      Logical NR,Micro_iter
      Character*80 Fullname
      
      Call Guesst2(T2old_aa,T2old_bb,T2old_ab,W0_aa,W0_bb,W0_ab,
     +             Fockoo_a,Fockoo_b,Fockvv_a,Fockvv_b,
     +             Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call Energy_frmt2(T2old_aa,T2old_bb,T2old_ab,W0_aa,W0_bb,W0_ab,
     +                  Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,E2_aa,E2_bb,E2_ab)
 
      Embpt2   = (E2_aa + E2_bb)*0.25D0 + E2_ab

      Write(6,"(a,F15.10)") "The MBPT(2) energy   : ", Embpt2

      Call Guesst1(T1old_aa,T1old_bb,Nocc_a,Nocc_b,Nvrt_a,
     +             Nvrt_b)

C Lets do a fractional non-HF MBPT(2) Brueckner calculations to
C define a set of orbitals for fractional orbital CC calculations 

      If (Brueck) Then
         Write(6,"(2a)") "Non-HF MBPT(2) Brueckner reference orbitals", 
     +                   " are obtained"
         Call Mbpt2_brueckner(T2old_aa,T2old_bb,T2old_ab,T1old_aa,
     +                        T1old_bb,W_aa,W_bb,W_ab,W0_aa,W0_bb,W0_ab,
     +                        W4_aa,W4_bb,W4_ab,W4_ba,W5_aa,W5_bb,W5_ab,
     +                        W5_ba,Fockoo_a,Fockoo_b,Fockvv_a,
     +                        Fockvv_b,Fockov_a,Fockov_b,Fme_a,Fme_b,
     +                        Embpt2,Work,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +                        Maxcor)
         Return 
      Endif 
      
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

      It1resid_aa = 1 
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
      NR        = .False.
      Micro_iter= .False.

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

C For Tau and Taup intermediate

         Call Form_tau(T2old_aa,T2old_bb,T2old_ab,Tau_aa,Tau_bb,
     +                 Tau_ab,Tap_aa,Tap_bb,Tap_ab,T1old_aa,T1old_bb,
     +                 Nocc_a,Nocc_b,Nvrt_a,Nvrt_b) 

C Fmi(m,i) intermediate 

         Call Form_fmi(Tap_aa,Tap_bb,Tap_ab,T1old_aa,T1old_bb,Fmi_a,
     +                 Fmi_b,W_aa,W_bb,W_ab,W5_aa,W5_bb,W5_ab,
     +                 W5_ba,Fockoo_od_a,Fockoo_od_b,Fockov_a,
     +                 Fockov_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
C Fae(a,e) intermediate 
      
         Call Form_fae(Tap_aa,Tap_bb,Tap_ab,T1old_aa,T1old_bb,Fae_a,
     +                 Fae_b,W_aa,W_bb,W_ab,W4_aa,W4_bb,W4_ab,
     +                 W4_ba,Fockvv_od_a,Fockvv_od_b,Fockov_a,
     +                 Fockov_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Fme(m,e) intermediate 

         Call Form_fme(Fme_a,Fme_b,T1old_aa,T1old_bb,W_aa,W_bb,
     +                 W_ab,Fockov_a,Fockov_b,Nocc_a,Nocc_b,
     +                 Nvrt_a,Nvrt_b,Fme_on)

C Form modified fae and Fmi intermediates (increment with Fme)
    
         Call form_modf_f(Fae_a,Fae_b,Fmi_a,Fmi_b,Fme_a,Fme_b,
     +                    T1old_aa,T1old_bb,Nocc_a,Nocc_b,
     +                    Nvrt_a,Nvrt_b,1.0D0) 

C Wmnij(mn,ij) intermediate 

         Call form_mnij(Tau_aa,Tau_bb,Tau_ab,T1old_aa,T1old_bb,
     +                  W1_aa,W1_bb,W1_ab,W_aa,W_bb,W_ab,W5_aa,
     +                  W5_bb,W5_ab,W5_ba,Wmnij_1,Wmnij_2,Wmnij_3,
     +                  Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,.False.)

C Wabef(ab,ef) intermediate 

         Call form_abef(Tau_aa,Tau_bb,Tau_ab,T1old_aa,T1old_bb,
     +                  W2_aa,W2_bb,W2_ab,W_aa,W_bb,W_ab,W4_aa,
     +                  W4_bb,W4_ab,W4_ba,Wabef_1,Wabef_2,Wabef_3,
     +                  Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,.False.)

C Wmbej(mb,ej) intermediate 

         Call form_mbej_t2(T2old_aa,T2old_bb,T2old_ab,W3_aaaa,
     +                     W3_bbbb,W3_abab,W3_baba,W3_abba,W3_baab,
     +                     W_aa,W_bb,W_ab,Wmbej_1,Wmbej_2,
     +                     Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,Nocc_a,
     +                     Nocc_b,Nvrt_a,Nvrt_b)

         Call form_mbej_t1(T1old_aa,T1old_bb,W_aa,W_bb,W_ab,W4_aa,
     +                     W4_bb,W4_ab,W4_ba,W5_aa,W5_bb,W5_ab,W5_ba,
     +                     Wmbej_1,Wmbej_2,Wmbej_3,Wmbej_4,Wmbej_5,
     +                     Wmbej_6,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Initialize the to <ab||ij>  and f(a,i)
 
         Call init_2abij(Work(It2resid_aa),Work(It2resid_bb),
     +                   Work(It2resid_ab),Work(It1resid_aa),
     +                   Work(It1resid_bb),W0_aa,W0_bb,W0_ab,
     +                   Fockov_a,Fockov_b,Nocc_a,Nocc_b,
     +                   Nvrt_a,Nvrt_b)

C Fae->T2

         Call Fae_int2(T2old_aa,T2old_bb,T2old_ab,Work(It2resid_aa),
     +                 Work(It2resid_bb),Work(It2resid_ab),Fae_a,
     +                 Fae_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Fmi->T2

         Call Fmi_int2(T2old_aa,T2old_bb,T2old_ab,Work(It2resid_aa),
     +                 Work(It2resid_bb),Work(It2resid_ab),Fmi_a,
     +                 Fmi_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C T1-> T2
         Call T1_int2(T1old_aa,T1old_bb,Work(It2resid_aa),
     +                Work(It2resid_bb),Work(It2resid_ab),W4_aa,
     +                W4_bb,W4_ab,W4_ba,W5_aa,W5_bb,W5_ab,W5_ba,
     +                Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Wmnij->T2
 
         Call Wmnij_int2(Tau_aa,Tau_bb,Tau_ab,Work(It2resid_aa),
     +                   Work(It2resid_bb),Work(It2resid_ab),Wmnij_1,
     +                   Wmnij_2,Wmnij_3,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Wabef->T2

         Call Wabef_int2(Tau_aa,Tau_bb,Tau_ab,Work(It2resid_aa),
     +                   Work(It2resid_bb),Work(It2resid_ab),Wabef_1,
     +                   Wabef_2,Wabef_3,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C T1T1->T2

         Call T1t1_int2(T1old_aa,T1old_bb,Work(It2resid_aa),
     +                  Work(It2resid_bb),Work(It2resid_ab),
     +                  W3_aaaa,W3_bbbb,W3_abab,W3_baba,W3_abba,
     +                  W3_baab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
         
C Wmbej-> T2

         Call Wmbej_int2(T2old_aa,T2old_bb,T2old_ab,Work(It2resid_aa),
     +                  Work(It2resid_bb),Work(It2resid_ab),Wmbej_1,
     +                  Wmbej_2,Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,
     +                  Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Form modified fae and Fmi intermediates (decrement with Fme)

         Call form_modf_f(Fae_a,Fae_b,Fmi_a,Fmi_b,Fme_a,Fme_b,
     +                    T1old_aa,T1old_bb,Nocc_a,Nocc_b,
     +                    Nvrt_a,Nvrt_b,-1.0D0)

C Fae->T1

         Call Fae_int1(T1old_aa,T1old_bb,Work(It1resid_aa),
     +                 Work(It1resid_bb),Fae_a,Fae_b,Nocc_a,Nocc_b,
     +                 Nvrt_a,Nvrt_b)

C Fmi->T1 
         Call Fmi_int1(T1old_aa,T1old_bb,Work(It1resid_aa),
     +                 Work(It1resid_bb),Fmi_a,Fmi_b,Nocc_a,Nocc_b,
     +                 Nvrt_a,Nvrt_b)

C Fme->T1
         Call Fme_int1(T2old_aa,T2old_bb,T2old_ab,Work(It1resid_aa),
     +                 Work(It1resid_bb),Fme_a,Fme_b,Nocc_a,Nocc_b,
     +                 Nvrt_a,Nvrt_b)

C T1->T1
         Call T1_int1(T1old_aa,T1old_bb,Work(It1resid_aa),
     +                Work(It1resid_bb),W3_aaaa,W3_bbbb,W3_abab,
     +                W3_baba,W3_abba,W3_baab,Nocc_a,Nocc_b,Nvrt_a,
     +                Nvrt_b)

C T2->T1 
         Call T2_int1(T2old_aa,T2old_bb,T2old_ab,W4_aa,W4_bb,W4_ab,
     +                W4_ba,W5_aa,W5_bb,W5_ab,W5_ba,Work(It1resid_aa),
     +                Work(It1resid_bb),Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
         If (NR) Then

         Call Fdiags_2r12(T1old_aa,T1old_bb,T2old_aa,T2old_bb,
     +                    T2old_ab,Work(It1resid_aa),
     +                    Work(It1resid_bb),Work(It2resid_aa),
     +                    Work(It2resid_bb),Work(It2resid_ab),
     +                    Fockoo_a,Fockoo_b,Fockvv_a,Fockvv_b,
     +                    Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

         Call Get_resid_stats(Work(It1resid_aa),Work(It1resid_bb),
     +                        Work(It2resid_aa),Work(It2resid_bb),
     +                        Work(It2resid_bb),Tresid_rms,
     +                        Tresid_max,Nocc_a,Nocc_b,Nvrt_a,
     +                        Nvrt_b) 

         Call form_modf_f(Fae_a,Fae_b,Fmi_a,Fmi_b,Fme_a,Fme_b,
     +                    T1old_aa,T1old_bb,Nocc_a,Nocc_b,
     +                    Nvrt_a,Nvrt_b,1.0D0)

         Call form_hbar(Fmi_a,Fmi_b,Fae_a,Fae_b,Fme_a,Fme_b,
     +                  Wmnij_1,Wmnij_2,Wmnij_3,
     +                  Wabef_1,Wabef_2,Wabef_3,Wmbej_1,Wmbej_2,
     +                  Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,
     +                  Wabci_1,Wabci_2,Wabci_3,Wabci_4,
     +                  Waibc_1,Waibc_2,Waibc_3,Waibc_4,
     +                  Wijka_1,Wijka_2,Wijka_3,Wijka_4,
     +                  Wiajk_1,Wiajk_2,Wiajk_3,Wiajk_4,
     +                  Wijab_1,Wijab_2,Wijab_3,
     +                  T1old_aa,T1old_bb,T2old_aa,T2old_bb,
     +                  T2old_ab,Tau_aa,Tau_bb,Tau_ab,W_aa,W_bb,W_ab,
     +                  W2_aa,W2_bb,W2_ab,W3_aaaa,W3_bbbb,W3_abab,
     +                  W3_baba,W3_abba,W3_baab,W4_aa,W4_bb,W4_ab,
     +                  W4_ba,W5_aa,W5_bb,W5_ab,W5_ba,Fockoo_a,
     +                  Fockoo_b,Fockvv_a,Fockvv_b,Fockov_a,
     +                  Fockov_b,Work(Iend),Memleft,Nocc_a,Nocc_b,
     +                  Nvrt_a,Nvrt_b)

         If (Micro_iter) Then

         Ihbar1d_aa = Iend
         Ihbar1d_bb = Ihbar1d_aa + L_aa
         Ihbar2d_aa = Ihbar1d_bb + L_bb
         Ihbar2d_bb = Ihbar2d_aa + L_aaaa
         Ihbar2d_ab = Ihbar2d_bb + L_bbbb
         Ismall_mat = Ihbar2d_ab + L_abab
         Iend       = Ismall_mat + Lineq_mxcyc*Lineq_mxcyc

         If (Iend .GT. Maxcor) Call Insmem("ccsd_light_driver",Iend,
     +                                      Maxcor)
         Memleft = Maxcor - Iend

         Call Hbar_diagonals(Fmi_a,Fmi_b,Fae_a,Fae_b,Wmnij_1,Wmnij_2,
     +                       Wmnij_3,Wabef_1,Wabef_2,Wabef_3,Wmbej_1,
     +                       Wmbej_2,Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,
     +                       Work(Ihbar1d_aa),Work(Ihbar1d_bb),
     +                       Work(Ihbar2d_aa),Work(Ihbar2d_bb),
     +                       Work(Ihbar2d_ab),Work(Iend),Memleft,
     +                       Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

         Call Get_new_r12(Fmi_a,Fmi_b,Fae_a,Fae_b,Fme_a,Fme_b,
     +                    Wmnij_1,Wmnij_2,Wmnij_3,
     +                    Wabef_1,Wabef_2,Wabef_3,Wmbej_1,Wmbej_2,
     +                    Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,
     +                    Wabci_1,Wabci_2,Wabci_3,Wabci_4,
     +                    Waibc_1,Waibc_2,Waibc_3,Waibc_4,
     +                    Wijka_1,Wijka_2,Wijka_3,Wijka_4,
     +                    Wiajk_1,Wiajk_2,Wiajk_3,Wiajk_4,
     +                    Wijab_1,Wijab_2,Wijab_3,Work(It1resid_aa),
     +                    Work(It1resid_bb),Work(It2resid_aa),
     +                    Work(It2resid_bb),Work(It2resid_ab),
     +                    Work(Ihbar1d_aa),Work(Ihbar1d_bb),
     +                    Work(Ihbar2d_aa),Work(Ihbar2d_bb),
     +                    Work(Ihbar2d_ab),Work(Ismall_mat),
     +                    T1old_aa,T1old_bb,T2old_aa,T2old_bb,
     +                    T2old_ab,Work(Iend),Memleft,Nocc_a,
     +                    Nocc_b,Nvrt_a,Nvrt_b)
         Else

         Ileng_ia_aa    = Nocc_a*Nvrt_a
         Ileng_ia_bb    = Nocc_b*Nvrt_b
         Ileng_ijab_aa  = Nvrt_a*Nvrt_a*Nocc_a*Nocc_a
         Ileng_ijab_bb  = Nvrt_b*Nvrt_b*Nocc_b*Nocc_b
         Ileng_ijab_ab  = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b

         Ileng_cc_jacobi = Ileng_ijab_aa+Ileng_ijab_bb+Ileng_ijab_ab+
     +                     Ileng_ia_aa+Ileng_ia_bb

         I_cc_jacobi = Iend
         Iend        = I_cc_jacobi + Ileng_cc_jacobi*Ileng_cc_jacobi
         Memleft     = Maxcor - Iend
         If (Iend .GT. Maxcor) Call Insmem("ccsd_light_driver",Iend,
     +                                      Maxcor)
     
         Call Form_jacobi(Fmi_a,Fmi_b,Fae_a,Fae_b,Fme_a,Fme_b,
     +                    Wmnij_1,Wmnij_2,Wmnij_3,
     +                    Wabef_1,Wabef_2,Wabef_3,Wmbej_1,Wmbej_2,
     +                    Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,
     +                    Wabci_1,Wabci_2,Wabci_3,Wabci_4,
     +                    Waibc_1,Waibc_2,Waibc_3,Waibc_4,
     +                    Wijka_1,Wijka_2,Wijka_3,Wijka_4,
     +                    Wiajk_1,Wiajk_2,Wiajk_3,Wiajk_4,
     +                    Wijab_1,Wijab_2,Wijab_3,Work(It1resid_aa),
     +                    Work(It1resid_bb),Work(It2resid_aa),
     +                    Work(It2resid_bb),Work(It2resid_ab),
     +                    T1old_aa,T1old_bb,T2old_aa,T2old_bb,
     +                    T2old_ab,Work(I_cc_jacobi),
     +                    Work(Iend),Memleft,Ileng_cc_jacobi,
     +                    Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
         Endif 

         Else

C T2old->T2new 

         Call Newt2(T2old_aa,T2old_bb,T2old_ab,Work(It2resid_aa),
     +              Work(It2resid_bb),Work(It2resid_ab),Fockoo_a,
     +              Fockoo_b,Fockvv_a,Fockvv_b,Nocc_a,Nocc_b,
     +              Nvrt_a,Nvrt_b)

C T1old->T1new 

         Call Newt1(T1old_aa,T1old_bb,Work(It1resid_aa),
     +              Work(It1resid_bb),Fockoo_a,Fockoo_b,Fockvv_a,
     +              Fockvv_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

         Endif 

         Call Energy_frmt2(T2old_aa,T2old_bb,T2old_ab,W0_aa,W0_bb,
     +                     W0_ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +                     E2ccsd_aa,E2ccsd_bb,E2ccsd_ab)

         Call Energy_frmt1(T1old_aa,T1old_bb,W0_aa,W0_bb,W0_ab,
     +                     Fockov_a,Fockov_b,Nocc_a,Nocc_b,Nvrt_a,
     +                     Nvrt_b,E1ccsd_aa,E1ccsd_bb,E1ccsd_ab,
     +                     Enhf_aa,Enhf_bb)

         Call Form_cumldens(T1old_aa,T1old_bb,Fockoo_a,Fockoo_b,
     +                      Work(Iend),Memleft,Nbasis,Ref_energy,
     +                      Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,Frac_occ)

         E2_new = (E2ccsd_aa + E2ccsd_bb)*0.25D0 + E2ccsd_ab
         E1_new = (E1ccsd_aa + E1ccsd_bb)*0.25D0 + E1ccsd_ab
         E3_new = Enhf_aa + Enhf_bb 
         E_new  = E2_new + E1_new + E3_new 
         E_diff= Dabs(E_new - E_old)
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
CSSS         If (Iter .Eq. 4) Converged = .True.

         IF (E_diff .LT. CC_conv) Converged = .True.
CSSS         IF (Tresid_rms .LT. CC_conv) Converged = .True.
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

#ifdef _DEBUG_LVL0
         Call Fdiags_2r12(T1old_aa,T1old_bb,T2old_aa,T2old_bb,
     +                    T2old_ab,Work(It1resid_aa),
     +                    Work(It1resid_bb),Work(It2resid_aa),
     +                    Work(It2resid_bb),Work(It2resid_ab),
     +                    Fockoo_a,Fockoo_b,Fockvv_a,Fockvv_b,
     +                    Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)


        Call form_modf_f(Fae_a,Fae_b,Fmi_a,Fmi_b,Fme_a,Fme_b,
     +                    T1old_aa,T1old_bb,Nocc_a,Nocc_b,
     +                    Nvrt_a,Nvrt_b,1.0D0)

        Call form_hbar(Fmi_a,Fmi_b,Fae_a,Fae_b,Fme_a,Fme_b,
     +                  Wmnij_1,Wmnij_2,Wmnij_3,
     +                  Wabef_1,Wabef_2,Wabef_3,Wmbej_1,Wmbej_2,
     +                  Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,
     +                  Wabci_1,Wabci_2,Wabci_3,Wabci_4,
     +                  Waibc_1,Waibc_2,Waibc_3,Waibc_4,
     +                  Wijka_1,Wijka_2,Wijka_3,Wijka_4,
     +                  Wiajk_1,Wiajk_2,Wiajk_3,Wiajk_4,
     +                  Wijab_1,Wijab_2,Wijab_3,
     +                  T1old_aa,T1old_bb,T2old_aa,T2old_bb,
     +                  T2old_ab,Tau_aa,Tau_bb,Tau_ab,W_aa,W_bb,W_ab,
     +                  W2_aa,W2_bb,W2_ab,W3_aaaa,W3_bbbb,W3_abab,
     +                  W3_baba,W3_abba,W3_baab,W4_aa,W4_bb,W4_ab,
     +                  W4_ba,W5_aa,W5_bb,W5_ab,W5_ba,Fockoo_a,
     +                  Fockoo_b,Fockvv_a,Fockvv_b,Fockov_a,
     +                  Fockov_b,Work(Iend),Memleft,Nocc_a,Nocc_b,
     +                  Nvrt_a,Nvrt_b)

         If (Micro_iter) Then

         Ihbar1d_aa = Iend
         Ihbar1d_bb = Ihbar1d_aa + L_aa
         Ihbar2d_aa = Ihbar1d_bb + L_bb
         Ihbar2d_bb = Ihbar2d_aa + L_aaaa
         Ihbar2d_ab = Ihbar2d_bb + L_bbbb
         Ismall_mat = Ihbar2d_ab + L_abab
         Iend       = Ismall_mat + Lineq_mxcyc*Lineq_mxcyc

         If (Iend .GT. Maxcor) Call Insmem("ccsd_light_driver",Iend,
     +                                      Maxcor)
         Memleft = Maxcor - Iend

         Call Hbar_diagonals(Fmi_a,Fmi_b,Fae_a,Fae_b,Wmnij_1,Wmnij_2,
     +                       Wmnij_3,Wabef_1,Wabef_2,Wabef_3,Wmbej_1,
     +                       Wmbej_2,Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,
     +                       Work(Ihbar1d_aa),Work(Ihbar1d_bb),
     +                       Work(Ihbar2d_aa),Work(Ihbar2d_bb),
     +                       Work(Ihbar2d_ab),Work(Iend),Memleft,
     +                       Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

         Call Get_new_r12(Fmi_a,Fmi_b,Fae_a,Fae_b,Fme_a,Fme_b,
     +                    Wmnij_1,Wmnij_2,Wmnij_3,
     +                    Wabef_1,Wabef_2,Wabef_3,Wmbej_1,Wmbej_2,
     +                    Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,
     +                    Wabci_1,Wabci_2,Wabci_3,Wabci_4,
     +                    Waibc_1,Waibc_2,Waibc_3,Waibc_4,
     +                    Wijka_1,Wijka_2,Wijka_3,Wijka_4,
     +                    Wiajk_1,Wiajk_2,Wiajk_3,Wiajk_4,
     +                    Wijab_1,Wijab_2,Wijab_3,Work(It1resid_aa),
     +                    Work(It1resid_bb),Work(It2resid_aa),
     +                    Work(It2resid_bb),Work(It2resid_ab),
     +                    Work(Ihbar1d_aa),Work(Ihbar1d_bb),
     +                    Work(Ihbar2d_aa),Work(Ihbar2d_bb),
     +                    Work(Ihbar2d_ab),Work(Ismall_mat),
     +                    T1old_aa,T1old_bb,T2old_aa,T2old_bb,
     +                    T2old_ab,Work(Iend),Memleft,
     +                    Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
         Else

         Ileng_ia_aa    = Nocc_a*Nvrt_a
         Ileng_ia_bb    = Nocc_b*Nvrt_b
         Ileng_ijab_aa  = Nvrt_a*(Nvrt_a-1)*Nocc_a*(Nocc_a-1)/4
         Ileng_ijab_bb  = Nvrt_b*(Nvrt_b-1)*Nocc_b*(Nocc_b-1)/4
         Ileng_ijab_ab  = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b

         Ileng_cc_jacobi = Ileng_ijab_aa+Ileng_ijab_bb+Ileng_ijab_ab+
     +                     Ileng_ia_aa+Ileng_ia_bb

         I_cc_jacobi = Iend
         Iend        = I_cc_jacobi + Ileng_cc_jacobi*Ileng_cc_jacobi
         Memleft     = Maxcor - Iend
         If (Iend .GT. Maxcor) Call Insmem("ccsd_light_driver",Iend,
     +                                      Maxcor)
         Call Form_jacobi(Fmi_a,Fmi_b,Fae_a,Fae_b,Fme_a,Fme_b,
     +                    Wmnij_1,Wmnij_2,Wmnij_3,
     +                    Wabef_1,Wabef_2,Wabef_3,Wmbej_1,Wmbej_2,
     +                    Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,
     +                    Wabci_1,Wabci_2,Wabci_3,Wabci_4,
     +                    Waibc_1,Waibc_2,Waibc_3,Waibc_4,
     +                    Wijka_1,Wijka_2,Wijka_3,Wijka_4,
     +                    Wiajk_1,Wiajk_2,Wiajk_3,Wiajk_4,
     +                    Wijab_1,Wijab_2,Wijab_3,Work(It1resid_aa),
     +                    Work(It1resid_bb),Work(It2resid_aa),
     +                    Work(It2resid_bb),Work(It2resid_ab),
     +                    T1old_aa,T1old_bb,T2old_aa,T2old_bb,
     +                    T2old_ab,Work(I_cc_jacobi),
     +                    Work(Iend),Memleft,Ileng_cc_jacobi,
     +                    Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

         Endif
#endif 

      Return
      End 

