      Subroutine Fix_hbar_abef(Hbar_abef1,Hbar_abef2,Hbar_abef3,
     +                         T2abef_1,T2_abef2,T2_abef3,
     +                         W_aa,W_bb,W_ab,T2aa,T2bb,T2ab,
     +                         Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z) 

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
 
      Dimension Hbar_abef1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_abef2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Hbar_abef3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension T2abef_1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension T2abef_2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension T2abef_3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Data Ione,One,Onem/1,1.0D0,-1.00D0/

      Lenabef_1 = Nvrt_a*Nvrt_a*Nvrt_a*Nvrt_a
      Lenabef_2 = Nvrt_b*Nvrt_b*Nvrt_b*Nvrt_b
      Lenabef_3 = Nvrt_a*Nvrt_b*Nvrt_a*Nvrt_b

C First form the 1/4t(ab,mn)<mn||cd> contribution.

      Call Form_t2abef(T2aa,T2bb,T2ab,W_aa,W_bb,W_ab,T2abef_1,
     +                 T2abef_2,T2abef_3,Nocc_a,Nocc_b,Nvrt_a,
     +                 Nvrt_b)

C Now add that to the existing Hbar to get 1/2t(ab,mn)<mn||cd> contribution

      Call Daxpy(Lenabef_1,One,T2abef_1,1,Hbar_abef1,1)
      Call Daxpy(Lenabef_2,One,T2abef_2,1,Hbar_abef2,1)
      Call Daxpy(Lenabef_3,One,T2abef_3,1,Hbar_abef3,1)


#ifdef _DEBUG_LVL0
      call checksum("Hbar_abef1:",Hbar_abef1,Lenabef_1)
      call checksum("Hbar_abef2:",Hbar_abef2,Lenabef_2)
      call checksum("Hbar_abef3:",Hbar_abef3,Lenabef_3)
#endif
      Return
      End
      
  
