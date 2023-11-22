      Subroutine Fix_hbar_ijkl(Hbar_ijkl1,Hbar_ijkl2,Hbar_ijkl3,
     +                         T2ijkl_1,T2_ijkl2,T2_ijkl3,
     +                         W_aa,W_bb,W_ab,T2aa,T2bb,T2ab,
     +                         Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z) 

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
 
      Dimension Hbar_ijkl1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Hbar_ijkl2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Hbar_ijkl3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension T2ijkl_1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension T2ijkl_2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension T2ijkl_3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Data Ione,One,Onem /1,1.0D0,-1.00D0/

C First form the 1/4t(ef,ij)<kl||ef> contribution.

      Lenijkl_1 = Nocc_a*Nocc_a*Nocc_a*Nocc_a
      Lenijkl_2 = Nocc_b*Nocc_b*Nocc_b*Nocc_b
      Lenijkl_3 = Nocc_a*Nocc_b*Nocc_a*Nocc_b

      Call Form_t2ijkl(T2aa,T2bb,T2ab,W_aa,W_bb,W_ab,T2ijkl_1,
     +                 T2ijkl_2,T2ijkl_3,Nocc_a,Nocc_b,Nvrt_a,
     +                 Nvrt_b)

C Now add that to the existing Hbar to get 1/2(ef,ij)<kl||ef> 

      Call Daxpy(Lenijkl_1,One,T2ijkl_1,1,Hbar_ijkl1,1)
      Call Daxpy(Lenijkl_2,One,T2ijkl_2,1,Hbar_ijkl2,1)
      Call Daxpy(Lenijkl_3,One,T2ijkl_3,1,Hbar_ijkl3,1)


#ifdef _DEBUG_LVL0
      call checksum("Hbar_ijkl1:",Hbar_ijkl1,Lenijkl_1)
      call checksum("Hbar_ijkl2:",Hbar_ijkl2,Lenijkl_2)
      call checksum("Hbar_ijkl3:",Hbar_ijkl3,Lenijkl_3)
#endif
      Return
      End
      
  
