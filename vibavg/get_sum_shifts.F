      Subroutine Get_sum_shifts(Delta,Equil,Nreals,Ndata) 

      Implicit Double Precision(A-H,O-Z)

      Dimension Delta(Nreals,Ndata)
      Dimension Equil(Nreals)

      Do K = 1, Ndata 
         Do J = 1, Nreals 
            Delta(J,K) = Delta(J,K) + Equil(J)
         Enddo
      Enddo
    
      Return
      End 
