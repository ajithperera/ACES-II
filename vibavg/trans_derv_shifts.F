      Subroutine Trans_derv_shifts(Prop,Factor,D,Nreals)

      Implicit Double Precision(A-H,O-Z)

      Dimension Prop(Nreals)
  
      Do I = 1, Nreals
         Prop(I) = (Factor*Prop(I))/D
      Enddo

      Return
      End 
