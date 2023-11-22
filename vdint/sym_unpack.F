      Subroutine Sym_unpack(Wfull, Wpack, Ndim, Nirrep, Nbfirr)

      Implicit Double Precision (A-H, O-Z)

      Dimension Wfull(Ndim, Ndim), Wpack(Ndim*(Ndim+1)/2), 
     &          Nbfirr(8)

      Idex = 0
      Kden = 0
     
      Do Irrep = 1, Nirrep

        Do I = 1 + Idex, Nbfirr(Irrep) + idex
           Do J = 1+idex, I
              Kden = Kden + 1
              If (I .Eq. J) Then
                 Wfull(j,i) = Wpack(Kden)
                 Wfull(i,j) = Wpack(Kden)
              Else
                 Wfull(j,i) = 0.50D0*Wpack(Kden)
                 Wfull(i,j) = 0.50D0*Wpack(Kden)
              Endif
         
           Enddo
        Enddo
           Idex = Idex + Nbfirr(Irrep)
      Enddo

      Return
      End
