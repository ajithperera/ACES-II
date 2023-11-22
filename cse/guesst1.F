      Subroutine Guesst1(T1aa_old,T1bb_old,Dens_aa,Dens_bb,Nocc_a,
     +                   Nocc_b,Nvrt_a,Nvrt_b,Nbasis)

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T1aa_old(Nvrt_a,Nocc_a)
      Dimension T1bb_old(Nvrt_b,Nocc_b)
      Dimension Dens_aa(Nbasis,Nbasis)
      Dimension Dens_bb(Nbasis,Nbasis)

      Integer I,A

#include "cse_vars.com"

#ifdef _NOSKIP
C AA block 

      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         If (A .EQ. I) Then
            T1aa_old(A,I) = Dens_aa(A,I) - Ocn_oa(i)
         Else
            T1aa_old(A,I) = Dens_aa(A,I) 
         Endif 
      ENDDO
      ENDDO

C BB block 

      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         If (A .EQ. I) Then
            T1bb_old(a,i) = Dens_bb(a,i) - Ocn_ob(i)
         Else 
            T1bb_old(a,i) = Dens_bb(a,i) 
         Endif 
      ENDDO
      ENDDO
#endif 

C#ifdef _NOSKIP
C AA block

      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         T1aa_old(A,I) = 0.0D0
      ENDDO
      ENDDO

C BB block

      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         T1bb_old(a,i) = 0.0D0
      ENDDO
      ENDDO
C#endif 

#ifdef _DEBUG_LVLM
      call output(Dens_aa,1,Nocc_a,1,Nvrt_a,Nocc_a,Nvrt_a,1) 
      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      call checksum("T1aa      :",T1aa_old,L_aa)
      call checksum("T1bb      :",T1bb_old,L_bb)
      do i=1,Nocc_a
      do a=1,Nvrt_a
      Write(6,"(2I2,1x,F15.10)") a,i,T1aa_old(a,i)
      enddo
      enddo
#endif 

      Return
      End
