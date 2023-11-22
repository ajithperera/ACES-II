      Subroutine Form_vv(T2cc,T2mp,T2,Iwork,Imaxcor,Nrow,Ncol,Iend,
     +                   Jend)

      Implicit Double Precision(A-H,O-Z)

      Dimension Iwork(Imaxcor)
      Dimension T2cc(Nrow,Ncol) 
      Dimension T2mp(Nrow,Ncol) 
      Dimension T2(Nrow,Ncol) 

      Do J = 1, Ncol
         Do I = 1, Nrow
            If (I .Le. Iend .and. J .Le. Jend) Then
               T2(I,J) = T2cc(I,J) 
            Endif
         Enddo
      Enddo 
   
      Return 
      End 
