      Subroutine Trans_derv_spncc(Prop,Factor,D,Nreals)

      Implicit Double Precision(A-H,O-Z)

      Dimension Prop(Nreals,Nreals) 
  
      Do J = 1, Nreals
         Do I = 1, Nreals
            Prop(I,J) = (Factor*Prop(I,J))/D
         Enddo
      Enddo 

      Return
      End 
