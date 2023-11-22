      Subroutine Form_rhs_t2rabij(T2aa_old,T2bb_old,T2ab_old,
     +                            Fmi_aa,Fmi_bb,Fae_aa,Fae_bb,
     +                            Wmnij_1,Wmnij_2,Wmnij_3,Wabef_1,
     +                            Wabef_2,Wabef_3,Wmbej_1,Wmbej_2,
     +                            Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,
     +                            T2resid_aa,T2resid_bb,T2resid_ab,
     +                            Fockoo_a,Fockoo_b,Fockvv_a,
     +                            Fockvv_b,Nocc_a,Nocc_b,Nvrt_a,
     +                            Nvrt_b)

      Implicit Double Precision (A-H,O-Z)
    
      Dimension Wmnij_1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Wmnij_2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Wmnij_3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension Wabef_1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Wabef_2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Wabef_3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension Wmbej_1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Wmbej_2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Wmbej_3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Wmbej_4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension Wmbej_5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension Wmbej_6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Dimension Fmi_aa(Nocc_a,Nocc_a)
      Dimension Fmi_bb(Nocc_b,Nocc_b)
      Dimension Fae_aa(Nvrt_a,Nvrt_a)
      Dimension Fae_bb(Nvrt_b,Nvrt_b)

      Dimension Fockoo_a(Nocc_a,Nocc_a)
      Dimension Fockoo_b(Nocc_b,Nocc_b)
      Dimension Fockvv_a(Nvrt_a,Nvrt_a)
      Dimension Fockvv_b(Nvrt_b,Nvrt_b)

      Dimension T2aa_old(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb_old(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab_old(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T2resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

#include "cse12_vars.com"

C Fae -> RHS_R2(ab,ij)

      Call T2fae_inrabij(T2aa_old,T2bb_old,T2ab_old,T2resid_aa,
     +                  T2resid_bb,T2resid_ab,Fae_aa,Fae_bb,
     +                  Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Fmi -> RHS_R2(ab,ij)

      Call T2fmi_inrabij(T2aa_old,T2bb_old,T2ab_old,T2resid_aa,
     +                  T2resid_bb,T2resid_ab,Fmi_aa,Fmi_bb,
     +                  Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C T2mnij-> RHS_R2(ab,ij)

      Call T2mnij_inrabij(T2aa_old,T2bb_old,T2ab_old,T2resid_aa,
     +                    T2resid_bb,T2resid_ab,Wmnij_1,Wmnij_2,
     +                    Wmnij_3,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C T2abef-> RHS_R2(ab,ij)

      Call T2abef_inrabij(T2aa_old,T2bb_old,T2ab_old,T2resid_aa,
     +                    T2resid_bb,T2resid_ab,Wabef_1,Wabef_2,
     +                    Wabef_3,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C T2mbej-> RHS_R2(ab,ij)

      Call T2mbej_inrabij(T2aa_old,T2bb_old,T2ab_old,T2resid_aa,
     +                    T2resid_bb,T2resid_ab,Wmbej_1,Wmbej_2,
     +                    Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,
     +                    Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

#ifdef _BUILT_RESIDUAL
      Call T2abij_inrabij(T2aa_old,T2bb_old,T2ab_old,T2resid_aa,
     +                    T2resid_bb,T2resid_ab,Fockoo_a,Fockoo_b,
     +                    Fockvv_a,Fockvv_b,Nocc_a,Nocc_b,Nvrt_a, 
     +                    Nvrt_b)
#endif 
      Return
      End 



