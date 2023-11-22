      Subroutine Scrnc_reord_cols(G2,G3,Nbfns,Irrep_ao,Iord)

      Implicit Double Precision(A-H,O-Z)

      Dimension G2(Nbfns*Nbfns,Nbfns*Nbfns)
      Dimension G3(Nbfns*Nbfns,Nbfns*Nbfns)
      Dimension Irrep_ao(Nbfns*Nbfns),Iordr_ao(Nbfns*Nbfns)

#include "syminf.com"
#include "sympop.com"
#include "sym.com"

      Nbfns2 = Nbfns * Nbfns 

      Ithru = 0
      Do Irrepr = 1, Nirrep
         Do I = 1, Nbfns2 
            If (Irrep_ao(I) .eq. Irrepr) Then
                Ithru = Ithru + 1
                Iordr_ao(Ithru) = I
            Endif
         Enddo
      Enddo

#ifdef _DEBUG_LVL0
      write(6,"(a,10(1x,I4))")"Reorder array      :", (Iordr_ao(i),
     +                       i=1,Nbfns2)
#endif
      Do I = 1, Nbfns2 
         Iold = Iordr_ao(I)
         Inew = I
            Do J = 1, Nbfns2
               G3(J,Inew) = G2(J,Iold)
            Enddo
      Enddo 
      Call Dcopy(Nbfns2*Nbfns2,G3,1,G2,1)

      Do I = 1, Nbfns2-1
         Do J = I+1, Nbfns2 
            If (Irrep_ao(I) .Gt. Irrep_ao(J)) Then
               Ikeep = Irrep_ao(J)
               Irrep_ao(J) = Irrep_ao(I)
               Irrep_ao(I) = Ikeep
#ifdef _NOSKIP
               Do K = 1, Nbfns2 
                  Dkeep = G2(K,J)
                  G2(K,J) = G2(K,I)
                  G2(K,I) = Dkeep
               Enddo
#endif 
            Endif 
         Enddo 
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The I<bra|ket>"
      Call output(G2,1,Nbfns2,1,Nbfns2,Nbfns2,Nnbfns2,1)
#endif 
    
      Return
      End
