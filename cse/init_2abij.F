      Subroutine Init_2abij(T1resid_aa,T1resid_bb,Fockov_aa,
     +                      Fockov_bb,Nocc_a,Nocc_b,Nvrt_a,
     +                      Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension T1resid_aa(Nvrt_a,Nocc_a)
      Dimension T1resid_bb(Nvrt_b,Nocc_b)
      Dimension Fockov_aa(Nocc_a,Nvrt_a) 
      Dimension Fockov_bb(Nocc_b,Nvrt_b) 

      Integer I,J,A,B

#include "cse_vars.com"

C AA block 

      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         C = Ocn_oa(I)*(1.0D0-Ocn_va(A))
         T1resid_aa(A,I) = fockov_aa(I,A)*C
#ifdef _DEBUG_LVLM
         Write(6,"(a,2I2,1x,3(1x,F15.10))")"A,I,T1resid,C1,C2: ",
     +                                      A,I,fockov_aa(A,I),
     +                                      Ocn_oa(I),(1.0D0-Ocn_va(A))
#endif 
      ENDDO
      ENDDO

C BB block

      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         C = Ocn_ob(i)*(1.0D0-Ocn_vb(a))
         T1resid_bb(a,i) = fockov_bb(i,a)*C 
      ENDDO
      ENDDO

#ifdef _DEBUG_LVL0
      Write(6,*) 
      L_aa   = Nocc_a*Nvrt_a
      L_bb   = Nocc_b*Nvrt_b
      call checksum("T1Resid_aa:",T1Resid_aa,L_aa)
      call checksum("T1Resid_bb:",T1Resid_bb,L_bb)
      do i=1,Nocc_a
      do a=1,Nvrt_a
      Write(6,"(2I2,1x,F15.10)") a,i,T1resid_aa(a,i)
      enddo
      enddo
#endif

      Return
      End
