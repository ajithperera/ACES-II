      Subroutine Print_cfc(Cfc,Nmodes)

      Implicit Double Precision(A-H,O-Z)

      Dimension Cfc(Nmodes,Nmodes,Nmodes)

      Do i=1,Nmodes
         Call output(Cfc(1,1,i),1,Nmodes,1,Nmodes,Nmodes,Nmodes,1)
      Enddo

      Return
      End
