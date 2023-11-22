      Subroutine Reord_eigs(Eigs,Ndim)

      implicit Double precision(A-H,O-Z)

      Dimension Eigs(Ndim)

      Do I = 1, Ndim-1
      Do J = I+1, Ndim
         If (Eigs(I) .Gt. Eigs(j)) Then
            Tmp = Eigs(I)
            Eigs(I) = Eigs(J)
            Eigs(J) = Tmp
         Endif 
      Enddo
      Enddo 
      
      Return
      End
