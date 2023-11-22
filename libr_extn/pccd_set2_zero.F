      Subroutine Pccd_set2_zero(T2,Nrow,Ncol)

      Implicit Double Precision(A-H,O-Z)

      Dimension T2(Nrow,Ncol)

      Data Dnull /0.0D0/

      DO K = 1, Ncol
         Do L = 1, Nrow
            If (K .Ne. L) T2(L,K) = DNull
         Enddo
      Enddo
      Return
      End
