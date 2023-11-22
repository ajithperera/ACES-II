      Subroutine Pccd_htau_trnspose(Htau_pq,Tmp,Nbas,Nocc,Nvrt)

      Implicit Double Precision (A-H,O-Z)

      Dimension Htau_pq(Nbas,Nbas) 
      Dimension Tmp(Nbas,Nbas) 

      Call Dcopy(Nbas*Nbas,Htau_pq,1,Tmp,1)

      Do Icol = Nocc+1, Nbas
         Do Irow = 1, Nocc
            Htau_pq(Icol,Irow) = Tmp(Irow,Icol)
         Enddo
      Enddo 

      Do Icol = Nocc+1, Nbas
         Do Irow = 1, Nocc
            Htau_pq(Irow,Icol) = Tmp(Icol,Irow)
         Enddo
      Enddo 

      Return
      End 
