      Subroutine Rccd_light_driver(W0_aa,W0_bb,W0_ab,W1_aa,W1_bb,
     +                             W1_ab,W2_aa,W2_bb,W2_ab,
     +                             W3_aaaa,W3_bbbb,W3_abab,W3_baba,
     +                             W3_abba,W3_baab,W_aa,W_bb,W_ab,
     +                             Wmnij_1,Wmnij_2,Wmnij_3,Wabef_1,
     +                             Wabef_2,Wabef_3,Wmbej_1,Wmbej_2,
     +                             Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,
     +                             T2old_aa,T2old_bb,T2old_ab,
     +                             Fockoo_a,Fockoo_b,Fockvv_a,
     +                             Fockvv_b,Fock_ov_a,Fock_ov_b,
     +                             Fockoo_od_a,Fockoo_od_b,
     +                             Fockvv_od_a,Fockvv_od_b,Fmi_a,
     +                             Fmi_b,Fae_a,Fae_b,Fme_a,Fme_b,
     +                             Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +                             Work,Maxcor)

      Implicit Double Precision(A-H,O-Z)

#include "ccdlight_vars.com"

      Dimension Work(Maxcor)

      Dimension T2old_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2old_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2old_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W0_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension W0_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
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

      Logical Converged,Fme_on

      Call Guesst2(T2old_aa,T2old_bb,T2old_ab,W0_aa,W0_bb,W0_ab,
     +             Fockoo_a,Fockoo_b,Fockvv_a,Fockvv_b,
     +             Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call Energy(T2old_aa,T2old_bb,T2old_ab,W0_aa,W0_bb,W0_ab,
     +            Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,E2_aa,E2_bb,E2_ab)
 
      Embpt2   = (E2_aa + E2_bb)*0.25D0 + E2_ab

      Write(6,"(a,F15.10)") "The MBPT(2) energy: ", Embpt2

      Call Symmetrize_t2(T2old_aa,T2old_bb,T2old_ab,Nocc_a,
     +                   Nocc_b,Nvrt_a,Nvrt_b)

      Ndim_oo_aa = Nocc_a * Nocc_a
      Ndim_oo_bb = Nocc_b * Nocc_b
      Ndim_oo_ab = Nocc_a * Nocc_b
      Ndim_vv_aa = Nvrt_a * Nvrt_a
      Ndim_vv_bb = Nvrt_b * Nvrt_b
      Ndim_vv_ab = Nvrt_a * Nvrt_b
      Ndim_ov_aa = Nocc_a * Nvrt_a
      Ndim_ov_bb = Nocc_b * Nvrt_b
      Ndim_ov_ab = Nocc_a * Nvrt_b
      Ndim_ov_ba = Nocc_b * Nvrt_a

      Iresid_aa = 1 
      Iresid_bb = Iresid_aa + Ndim_oo_aa*Ndim_vv_aa
      Iresid_ab = Iresid_bb + Ndim_oo_bb*Ndim_vv_bb
      Iend      = Iresid_ab + Ndim_oo_ab*Ndim_vv_ab
  
      If (Iend .GT. Maxcor) Call Insmem("rccd_light_driver",Iend,
     +                                   Maxcor)
      Iter      = 0
      Converged = .False. 
      E_old     = Embpt2
      E_corr(0) = Embpt2 

      Do while ((.NOT. Converged) .AND. Iter .Le. Cc_maxcyc) 
         Iter = Iter + 1

         Call form_mbej_rc(T2old_aa,T2old_bb,T2old_ab,W3_aaaa,
     +                     W3_bbbb,W3_abab,W3_baba,W3_abba,W3_baab,
     +                     W_aa,W_bb,W_ab,Wmbej_1,Wmbej_2,
     +                     Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,Nocc_a,
     +                     Nocc_b,Nvrt_a,Nvrt_b)

C Initialize the to <ab||ij> 

         Call init_2abij(Work(Iresid_aa),Work(Iresid_bb),
     +                   Work(Iresid_ab),W0_aa,W0_bb,W0_ab,
     +                   Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Wmbej-> T2

        Call Wmbej_cont_rc(T2old_aa,T2old_bb,T2old_ab,Work(Iresid_aa),
     +                     Work(Iresid_bb),Work(Iresid_ab),Wmbej_1,
     +                     Wmbej_2,Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,
     +                     Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C T2old->T2new 

         Call Newt2(T2old_aa,T2old_bb,T2old_ab,Work(Iresid_aa),
     +              Work(Iresid_bb),Work(Iresid_ab),Fockoo_a,
     +              Fockoo_b,Fockvv_a,Fockvv_b,Nocc_a,Nocc_b,
     +              Nvrt_a,Nvrt_b)

         Call Energy(T2old_aa,T2old_bb,T2old_ab,W0_aa,W0_bb,W0_ab,
     +               Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,Eccd_aa,Eccd_bb,
     +               Eccd_ab)

         E_new =(Eccd_aa + Eccd_bb)*0.25D0 + 0.50D0*Eccd_ab
         E_diff= Dabs(E_new - E_old)
         E_corr(Iter) = E_new 

#ifdef _DEBUG_LVL0
         Write(6,"(a,2(1x,F15.10),2x,a)")
     +         "The energy and the energy difference:", E_new, E_diff,
     +          "at"
         Write(6,"(a,1x,I3))") "iteration no.: ", Iter
#endif 
CSSS         If (Iter .Eq. 2) Converged = .True.
         IF (E_diff .LT. CC_conv) Converged = .True.
 
         E_old = E_new 

      enddo

      Call Getrec(20,"JOBARC","SCFENEG ",1,E_scf)

      If (Converged) Then  
         Write(6,*)
         Write(6,9) 
         Write(6,10)
         Write(6,11)
         Write(6,10)
      Else
         Write(6,*)
         Write(6,"(a,1x,I3,a)") " The Ring-CCD did not converged in", 
     +                          cc_maxcyc, " iterations."
      Endif 

 9    Format(T10,'Summary of iterative solution of Ring-CCD equations')
 10   Format(6x,59("-")) 
 11   Format(T24,'Correlation',T55,'Total',/,T8,'Iteration',T26,
     +       'Energy',T55,'Energy')
      Do I = 0, Iter
         Write(6,"(5x,I3,10x,F18.12,10x,F19.12)") I, E_corr(i), 
     +                                             E_corr(i)+E_scf
      Enddo
      Write(6,10)
C           -          -
C Construct F(a,e) and F(m,i) to check whether I can match them in the
C general code.

      Call Form_fmi_rc(T2old_aa,T2old_bb,T2old_ab,Fmi_a,Fmi_b,
     +                 W_aa,W_bb,W_ab,Fockoo_od_a,Fockoo_od_b,
     +                 Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call Form_fae_rc(T2old_aa,T2old_bb,T2old_ab,Fae_a,Fae_b,
     +                 W_aa,W_bb,W_ab,Fockvv_od_a,Fockvv_od_b,
     +                 Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Ndim_oo_aa = Nocc_a * Nocc_a
      Ndim_oo_bb = Nocc_b * Nocc_b
      Ndim_oo_ab = Nocc_a * Nocc_b
      Ndim_vv_aa = Nvrt_a * Nvrt_a
      Ndim_vv_bb = Nvrt_b * Nvrt_b
      Ndim_vv_ab = Nvrt_a * Nvrt_b
      Ndim_ov_aa = Nocc_a * Nvrt_a
      Ndim_ov_bb = Nocc_b * Nvrt_b
      Ndim_ov_ab = Nocc_a * Nvrt_b
      Ndim_ov_ba = Nocc_b * Nvrt_a

      Ihmbej_aaaa = Iend 
      Ihmbej_bbbb = Ihmbej_aaaa + Ndim_ov_aa*Ndim_ov_aa
      Ihmbej_abab = Ihmbej_bbbb + Ndim_ov_bb*Ndim_ov_bb
      Ihmbej_baba = Ihmbej_abab + Ndim_ov_ab*Ndim_ov_ba
      Ihmbej_abba = Ihmbej_baba + Ndim_ov_ba*Ndim_ov_ab
      Ihmbej_baab = Ihmbej_abba + Ndim_oo_aa*Ndim_vv_bb
      Iend        = Ihmbej_baab + Ndim_oo_bb*Ndim_vv_aa

      If (Iend .GT. Maxcor) Call Insmem("rccd_light_driver",Iend,
     +                                   Maxcor)

      Call form_mbej_rc_hbar(T2old_aa,T2old_bb,T2old_ab,W3_aaaa,
     +                       W3_bbbb,W3_abab,W3_baba,W3_abba,W3_baab,
     +                       W_aa,W_bb,W_ab,Work(Ihmbej_aaaa),
     +                       Work(Ihmbej_bbbb),Work(Ihmbej_abab),
     +                       Work(Ihmbej_baba),Work(Ihmbej_abba),
     +                       Work(Ihmbej_baab),Nocc_a,Nocc_b,
     +                       Nvrt_a,Nvrt_b)
      
      Return
      End 

