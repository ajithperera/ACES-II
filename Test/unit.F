       Subroutine Unit(D,NMO)

       Implicit Double Precision (A-H,O-Z)

       Dimension D(NMO,NMO) 

        
       Do i=1, Nmo
       Do j=1, Nmo
          if (i .eq. j) then
            D(j,i) = 1.0D0
          Else
            D(j,i) = 0.0D0
          Endif 
       Enddo
       Enddo 

       Return
       End
