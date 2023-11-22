
      Subroutine test(
     +                Fock_ov_a,Fock_ov_b,
     +                Fockoo_od_a,Fockoo_od_b,
     +                Fockvv_od_a,Fockvv_od_b,
     +                Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +                Work,Maxcor)

      Implicit Double Precision (A-H,O-Z)

#include "ccsdlight_vars.com"

      Dimension Work(Maxcor)

      Dimension fockov_a(Nocc_a,Nvrt_a)
      Dimension fockov_b(Nocc_b,Nvrt_b)

      Dimension fockoo_od_a(Nocc_a,Nocc_a)
      Dimension fockoo_od_b(Nocc_b,Nocc_b)
      Dimension fockvv_od_a(Nvrt_a,Nvrt_a)
      Dimension fockvv_od_b(Nvrt_b,Nvrt_b)

      call check_fov(Fockov_a,Fockov_b,Nocc_a,Nocc_b,Nvrt_a,
     +                 Nvrt_b)

#ifdef _NOSKIP

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

        Call Hbar_t_mult(Fmi_a,Fmi_b,Fae_a,Fae_b,Fme_a,Fme_b,
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
     +                    T2old_ab,Work(Iend),Memleft,Nocc_a,Nocc_b,
     +                    Nvrt_a,Nvrt_b)

#endif




      Return
      End

