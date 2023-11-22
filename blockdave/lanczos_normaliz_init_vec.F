      Subroutine Lanczos_normaliz_init_vec(Data_in_out,Nsize)

      Implicit Double Precision (A-H,O-Z)
      
      Dimension Data_in_out(Nsize)
      Double Precision Norm 
      
      Norm = Ddot(Nsize,Data_in_out,1,Data_in_out,1)

      If (Norm .Ne. 0.0D0) Then
         Sqrt_Norm = 1.0D0/Dsqrt(Norm)
      Else
         Write(6,"(a)") " Division by zero @Tdee_normaliz_init_vec"
         call Errex 
      Endif 

      Call Dscal(Nsize,Sqrt_Norm,Data_in_out,1)
  
      Return
      End
