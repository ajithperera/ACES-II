      Subroutine Form(T1cc,T1,Iwork,Imaxcor,Nrow,Ncol,Aend,Iend)

      Implicit Double Precision(A-H,O-Z)

      Integer Aend 
      Dimension Iwork(Imaxcor)
      Dimension T1cc(Nrow,Ncol) 
      Dimension T1(Nrow,Ncol) 
     
      Jend = Ncol - (Iend-1)

      Do J = Ncol, Jend, -1 
         Do I = 1, Nrow 
            If (I .Le. Aend) Then
               T1(I,J) = T1cc(I,J) 
            Endif
         Enddo
      Enddo 
   
      Return 
      End 
