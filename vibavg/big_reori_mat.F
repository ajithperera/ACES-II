      Subroutine Big_reori_mat(Rbig,R,Nreals)

      Implicit Double Precision(A-H,O-Z)

      Dimension R(3,3), Rbig(3*Nreals,3*Nreals)


      Call Dzero(Rbig,9*Nreals*Nreals)

      Ioff1 = 0
      Ioff2 = 0

      Do Iatm = 1, Nreals
         Do Ixyz = 1, 3
            Do Jxyz = 1, 3
               RBig(Ioff1+Ixyz,Ioff2+Jxyz) = R(Ixyz,Jxyz)
            Enddo
         Enddo
         Ioff1 = Ioff1 + 3
         Ioff2 = Ioff2 + 3
      Enddo

      Return
      End 
  
