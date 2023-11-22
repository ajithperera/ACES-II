      Subroutine Form_rhs_t1t1rabij(T2aa_old,T2bb_old,T2ab_old,
     +                            Wmnij_1,Wmnij_2,Wmnij_3,Wabef_1,
     +                            Wabef_2,Wabef_3,T2resid_aa,
     +                            T2resid_bb,T2resid_ab,Nocc_a,
     +                            Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)
    
      Dimension Wmnij_1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Wmnij_2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Wmnij_3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension T2aa_old(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb_old(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab_old(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Wabef_1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Wabef_2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Wabef_3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension T2resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

#include "cse12_vars.com"

C T2mnij-> RHS_R2(ab,ij)

      Call T2mnij_inrabij(T2aa_old,T2bb_old,T2ab_old,T2resid_aa,
     +                    T2resid_bb,T2resid_ab,Wmnij_1,Wmnij_2,
     +                    Wmnij_3,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C T2abef-> RHS_R2(ab,ij)

      Call T2abef_inrabij(T2aa_old,T2bb_old,T2ab_old,T2resid_aa,
     +                    T2resid_bb,T2resid_ab,Wabef_1,Wabef_2,
     +                    Wabef_3,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)


      Return
      End 



