      Subroutine Form_residuals(T1aa,T1bb,T2aa,T2bb,T2ab,R1ai_aa,
     +                          R1ai_bb,R2abij_aa,R2abij_bb,
     +                          R2abij_ab,Fockoo_a,Fockoo_b,Fockvv_a,
     +                          Fockvv_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +                          Rms) 

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)
      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension R1ai_aa(Nvrt_a,Nocc_a)
      Dimension R1ai_bb(Nvrt_b,Nocc_b)
      Dimension R2abij_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension R2abij_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension R2abij_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Fockoo_a(Nocc_a,Nocc_a)
      Dimension Fockoo_b(Nocc_b,Nocc_b)
      Dimension Fockvv_a(Nvrt_a,Nvrt_a)
      Dimension Fockvv_b(Nvrt_b,Nvrt_b)

      Dimension V1_aa(6),V1_bb(6),V2_aa(6),V2_bb(6),V2_ab(6)

#include "cse12_vars.com"

      Do = 1.0D0
      Call Form_rai(T1aa,T1bb,R1ai_aa,R1ai_bb,Fockoo_a,Fockoo_b,
     +              Fockvv_a,Fockvv_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +              Do)

      Call Form_rabij(T2aa,T2bb,T2ab,R2abij_aa,R2abij_bb,R2abij_ab,
     +                Fockoo_a,Fockoo_b,Fockvv_a,Fockvv_b,Nocc_a,
     +                Nocc_b,Nvrt_a,Nvrt_b,Do)

      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b

      Call Vstat(R1ai_aa,V1_aa,L_aa)
      Call Vstat(R1ai_bb,V1_bb,L_bb)
      Call Vstat(R2abij_aa,V2_aa,L_aaaa)
      Call Vstat(R2abij_bb,V2_bb,L_bbbb)
      Call Vstat(R2abij_ab,V2_ab,L_abab)

      Rms = Max(V1_aa(5),V1_bb(5))
      Rms = Max(Rms,V2_aa(5))
      Rms = Max(Rms,V2_bb(5))
      Rms = Max(Rms,V2_ab(5))

#ifdef _DEBUG_LVL0
      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("R1ai_aa   :",R1ai_aa,L_aa)
      call checksum("R1ai_bb   :",R1ai_bb,L_bb)
      call checksum("R2abij_aa :",R2abij_aa,L_aaaa)
      call checksum("R2abij_bb :",R2abij_bb,L_bbbb)
      call checksum("R2abij_ab :",R2abij_ab,L_abab)
      Write(6,"(a,F20.12)") "The RMS              :", RMS
      Write(6,*) 
#endif
      Undo = -1.0D0
      Call Form_rai(T1aa,T1bb,R1ai_aa,R1ai_bb,Fockoo_a,Fockoo_b,
     +              Fockvv_a,Fockvv_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +              Undo)

      Call Form_rabij(T2aa,T2bb,T2ab,R2abij_aa,R2abij_bb,R2abij_ab,
     +                Fockoo_a,Fockoo_b,Fockvv_a,Fockvv_b,Nocc_a,
     +                Nocc_b,Nvrt_a,Nvrt_b,Undo)

      Return
      End
