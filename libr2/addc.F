      Subroutine addc(Length, Work, constant)

      Implicit double precision(A-H, O-Z)

      Dimension Work(Length)
    
      DO I =1, Length

       Work(i) = Work(i) + Constant
    
      Enddo
  
      Return 
      End
