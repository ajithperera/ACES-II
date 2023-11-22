      subroutine Set_active_space(Work,Maxcor,Nocc_a,Nocc_b,Nvrt_a,
     +                            Nvrt_b,Length)
     
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Integer Act_start

#include "ccsdlight_vars.com"

      L1_aa  = Nocc_a*Nvrt_a
      L1_bb  = Nocc_b*Nvrt_b
      L2_aa  = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L2_bb  = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L2_ab  = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      Length = L1_aa+L1_bb+L2_aa+L2_bb+L2_ab

      Ibgn  = 1
      Iend  = Ibgn + Length 
   
      If (Iend .Gt. Maxcor)
     +   Call Insmem("@-set_active_space",Iend,Maxcor)
     
      Call Form_exc_mask(Work(Ibgn),Length,Nocc_a,Nocc_b,Nvrt_a,
     +                   Nvrt_b)

#ifdef _DEBUG_LVL0
      Sum = 0.0D0
      Do i=1,Length
         Sum = Sum + Work(Ibgn+i-1)*Work(Ibgn+i-1)
      Enddo
      Write(6,"(a,1x,2(1x,ES8.2E2))") " @-set_active_space,P*P     = ",
     +                                  Sum
      Write(6,*) 
#endif

      Return
      End
