      Subroutine Pccd_sum(Hpq,M,E,S)

      Implicit Double Precision(A-H,O-Z)
      Character*1 S
      Logical Symmetry 
      Dimension Hpq(M,M)

      Common /Symm/Symmetry 
      call output(Hpq,1,M,1,M,M,M,1)
      If (S .EQ. "D") Then
      Do I = 1, M
         E = E + Hpq(I,I)
      Enddo
      Else IF (S .EQ. "F") THEN
      Do I = 1, M
      Do J = 1, M
         E = E + Hpq(I,J)
      Enddo
      Enddo
      Endif

      Return
      End
  
     
