      Subroutine Symmetrize_t2(T2aa_old,T2bb_old,T2ab_old,Nocc_a,
     +                         Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T2aa_old(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb_old(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab_old(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Integer I,J,A,B

#include "ccdlight_vars.com"

C AAAA and BBBB block 

      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
       
      Call Dcopy(L_aaaa,T2ab_old,1,T2aa_old,1)
      Call Dcopy(L_bbbb,T2ab_old,1,T2bb_old,1)

#ifdef _DEBUG_LVL0
      Write(6,*) 
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2aa      :",T2aa_old,L_aaaa)
      call checksum("T2bb      :",T2bb_old,L_bbbb)
      call checksum("T2ab      :",T2ab_old,L_abab)
#endif

      Return
      End
