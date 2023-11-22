      Subroutine Form_rhs_h2trabij(T2aa_old,T2bb_old,T2ab_old,
     +                             Htmi_aa,Htmi_bb,Htae_aa,Htae_bb,
     +                             Htmnij_1,Htmnij_2,Htmnij_3,Htabef_1,
     +                             Htabef_2,Htabef_3,Htmbej_1,Htmbej_2,
     +                             Htmbej_3,Htmbej_4,Htmbej_5,Htmbej_6,
     +                             Htabij_1,Htabij_2,Htabij_3,
     +                             T2resid_aa,T2resid_bb,T2resid_ab,
     +                             Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)
    
      Dimension Htmnij_1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Htmnij_2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Htmnij_3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension Htabef_1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Htabef_2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Htabef_3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension Htabij_1(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Htabij_2(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Htabij_3(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Htmbej_1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Htmbej_2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Htmbej_3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Htmbej_4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension Htmbej_5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension Htmbej_6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Dimension Htmi_aa(Nocc_a,Nocc_a)
      Dimension Htmi_bb(Nocc_b,Nocc_b)
      Dimension Htae_aa(Nvrt_a,Nvrt_a)
      Dimension Htae_bb(Nvrt_b,Nvrt_b)

      Dimension T2aa_old(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb_old(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab_old(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T2resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

#include "cse12_vars.com"

C Htmnij-> RHS_R2(ab,ij)

      Call Htmnij_inrabij(T2aa_old,T2bb_old,T2ab_old,T2resid_aa,
     +                    T2resid_bb,T2resid_ab,Htmnij_1,Htmnij_2,
     +                    Htmnij_3,Htabij_1,Htabij_2,Htabij_3,
     +                    Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
C Htabef-> RHS_R2(ab,ij)

      Call Htabef_inrabij(T2aa_old,T2bb_old,T2ab_old,T2resid_aa,
     +                     T2resid_bb,T2resid_ab,Htabef_1,Htabef_2,
     +                     Htabef_3,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Htae -> RHS_R2(ab,ij)

      Call Htae_inrabij(T2aa_old,T2bb_old,T2ab_old,T2resid_aa,
     +                  T2resid_bb,T2resid_ab,Htae_aa,Htae_bb,
     +                  Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Htmi -> RHS_R2(ab,ij)

      Call Htmi_inrabij(T2aa_old,T2bb_old,T2ab_old,T2resid_aa,
     +                  T2resid_bb,T2resid_ab,Htmi_aa,Htmi_bb,
     +                  Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Htmbej-> RHS_R2(ab,ij)

      Call Htmbej_inrabij(T2aa_old,T2bb_old,T2ab_old,T2resid_aa,
     +                    T2resid_bb,T2resid_ab,Htmbej_1,Htmbej_2,
     +                    Htmbej_3,Htmbej_4,Htmbej_5,Htmbej_6,
     +                    Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Return
      End 



