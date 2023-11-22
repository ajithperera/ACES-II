      Subroutine Sympack(Wfull, Wpack, Ndim, Nirrep, Nbfirr)

      Implicit Double Precision (A-H, O-Z)

      Dimension Wfull(Ndim, Ndim), Wpack(Ndim*(Ndim+1)/2), 
     &          Nbfirr(8)

      Idex = 0
      Kden = 0
      Do Irrep = 1, Nirrep

        Do I = 1 + Idex, Nbfirr(Irrep) + idex
           Do J = 1+idex, I
              Kden = Kden + 1
              Wpack(Kden) = Wfull(j,i)
           Enddo
        Enddo
           Idex = Idex + Nbfirr(Irrep)
      Enddo

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "@-Sympack, The symmetry packed array"
      Write(6,*) idex
      Write(6, "(6(1x,F10.7))") (Wpack(I), I=1, Idex*(Idex+1)/2)
#endif

      Return
      End
