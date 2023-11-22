      Subroutine Pccd_check_htau(Label,Hpq,N,Type,S)

      Implicit Double Precision(A-H,O-Z)
      Character*9 Label 
      Character*2 Type
      Character*1 S
      Logical Symmetry
      Dimension Hpq(N*N)

#include "syminf.com"   
#include "sym.com"   
      Common /Symm/Symmetry 

      Data Dnull,Ione,Inull /0.0D0,1,0/

      If (Type .Eq. "OO" .OR. Type .Eq. "VV") THEN
         E   = Dnull
         Ioff= Ione
         write(6,"(a,2(2x,F15.10))") Label
CSSS      write(6,"(a,2(2x,F15.10))") Label, e,e**2 
         do Irrep = 1, Nirrep
            If (Symmetry) Write(6,"(a,I2)") "Irrep:",Irrep
            If (Type .Eq. "OO") M = Pop(Irrep,1)
            If (Type .Eq. "VV") M = Vrt(Irrep,1)
            Call Pccd_sum(Hpq(Ioff),M,E,S)
            Ioff = Ioff + M*M
         enddo
      Else 
         Ioff= Ione
         do Irrep = 1, Nirrep
            If (Type .Eq. "OV") Then
                If (Symmetry) Write(6,"(a,I2)") "Irrep:",Irrep
                M = Pop(Irrep,1)
                N = Vrt(Irrep,1)
                Call output(Hpq(Ioff),1,M,1,N,M,N,1)
               Ioff= Ioff + M*N
            Elseif (Type .Eq. "VO") Then
                If (Symmetry) Write(6,"(a,I2)") "Irrep:",Irrep
	        M = Vrt(Irrep,1)
	        N = Pop(Irrep,1)
	        Call output(Hpq(Ioff),1,M,1,N,M,N,1)
                Ioff= Ioff + M*N
            Endif 
         enddo
      Endif 
 
      return
      end 
