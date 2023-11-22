      Subroutine Get_sum_spincc(Delta,Equil,Nreals,Ndata) 

      Implicit Double Precision(A-H,O-Z)

      Dimension Delta(Nreals,Nreals,Ndata)
      Dimension Equil(Nreals,Nreals)

      Do K = 1, Ndata 
         Do J = 1, Nreals 
            Do I = Nreals,J+1,-1
               Delta(I,J,K) = Delta(I,J,K) + Equil(I,J)
            Enddo
         Enddo
      Enddo
    
      Return
      End 
