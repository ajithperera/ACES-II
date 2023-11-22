      Subroutine Build_reori_matrix(Qnew,Qold,AtmMass,W,Natoms)

      Implicit Double Precision(A-H, O-Z)

#include "mxatms.par"
      Dimension Qnew(3,Natoms) 
      Dimension Qold(3,Natoms) 
      Dimension AtmMass(Natoms) 
      Dimension Scr(6*Mxatms+18)
      Dimension W(3,3)

      Data Done /1.0D0/

#ifdef _DEBUG_LVL0
      Print*, Natoms 
      Write(6,"(a)") " The new (sym. rotated) geometry"
      Do I = 1, Natoms
         Write(6,"(3(1x,F12.6))") (Qnew(j,i),j=1,3)
      Enddo 
      Write(6,"(a)") " The old (principal axis) geometry"
      Do I = 1, Natoms
         Write(6,"(3(1x,F12.6))") (Qold(j,i),j=1,3)
      Enddo 
#endif 
C Qucik return for atoms and diatomics. The transfornmation is 
C the unit matrix.

      If (Natoms .Eq. 1 .OR. Natoms .Eq. 2) Then
         Call Zero(W,9)
         W(1,1) = Done 
         W(2,2) = Done 
         W(3,3) = Done  
         Return 
      Endif 

      Call Q2qprime(Qnew,Qold,AtmMass,Scr,W,Natoms)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " Q(Sym. ordered)<->Q(Principle axis ordered) "
      Do j = 1,3
      Write(6,"(3(1x,F12.6))") (W(j,i),i=1,3)
      Enddo 
#endif 

      Return
      End
