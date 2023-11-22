      Subroutine Get_cfc_contr(Cfc,Cfc_contr,Exponent,Omega,Temp,
     +                         Ndata,Nmodes,Nreals,Icol)

      Implicit Double Precision(A-H,O-Z)
      Double Precision Num

#include "const.h"
      Dimension Cfc(Nmodes,Nmodes,Nmodes)
      Dimension Cfc_contr(Nmodes,Ndata)
      Dimension Omega(3*Nreals)
      Dimension Exponent(Nmodes,Ndata)
      Dimension Temp(Ndata)

      Data Done,Two,Tensq/1.0D0,2.0D0,100.0D0/
      
C Omega (frequncies in cm-1),C=in cm s^-1,h is in Js,K_b is in 
C J K^-1 and Temp. in K. 

#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,"(a)")" Printing from get_cfc_contr"
      write(6,"(a)") " The temperatures:"
      Write(6,"(4(1x,F10.4))") (Temp(i),i=1,Ndata)
      Write(6,"(a)") " The harmonic frequencies:"
      Write(6,"(4(1x,F15.6))") (Omega(i),i=1,3*Nreals)
#endif 

      Do J = 1, Ndata
         Do I = 1, Nmodes
            Num = (Omega(I+Icol)*CONST_h*CONST_c*Tensq)
            Den = (Two*Temp(J)*CONST_k)
            Exp = Num/Den
            Exponent(I,J) = Done/Dtanh(Exp)
         Enddo
      Enddo

C Cubic force constants are in cm-1 units and then property derivatives 
C are in their natural units. 

      Do J = 1, Ndata 
         Do I = 1, Nmodes 
            Do K = 1, Nmodes 
               Cfc_contr(I,J) = Cfc_contr(I,J) + Cfc(I,K,K)*
     +                          Exponent(K,J)
             Enddo 
         Enddo 
      Enddo 

#ifdef _DEBUG_LVL0
      Write(6,*) 
      write(6,"(a)") " The cubic force constant contribution:"
      do j=1,Ndata
      Write(6,"(6(1x,ES18.9E2))") (Cfc_contr(I,J),I=1,Nmodes)
      enddo
      write(6,"(a)") " The exponent contribution:"
      Do J=1, Ndata
      Write(6,"(6(1x,ES18.9E2))") (Exponent(I,J),I=1,Nmodes)
      Enddo 
      Write(6,*)
#endif

      Return
      End 
     
