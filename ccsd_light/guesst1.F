      Subroutine Guesst1(T1aa_old,T1bb_old,Nocc_a,Nocc_b,Nvrt_a,
     +                   Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T1aa_old(Nvrt_a,Nocc_a)
      Dimension T1bb_old(Nvrt_b,Nocc_b)

      Integer I,A

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

      Return
      End
