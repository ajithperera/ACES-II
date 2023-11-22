      Subroutine Scrnc_sympack_aoints(W_in,W_out,Nbfns,Irrepx)

      Implicit Double Precision(A-H,O-Z)
      Dimension W_in(Nbfns,Nbfns), W_out(Nbfns*Nbfns)

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)
#include "syminf.com"
#include "sympop.com"
#include "sym.com"

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The symmetry unpacked AO integrals"
      call output(W_in,1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
#endif 
      Ndim2 = Nbfns * Nbfns 
      Call Dzero(W_out,Ndim2)

      Ithru=0
      Do Irrepr=1, Nirrep
         Irrepl=Dirprd(Irrepr,Irrepx)
          Do Nu=1, Iaopop(Irrepr)
             Do Mu=1,Iaopop(Irrepl) 
                Ithru = Ithru + 1
                Mu_off = Mu + IoffAO(Irrepl)-1 
                Nu_off = Nu + IoffAO(Irrepr)-1 
                W_out(Ithru) = W_in(Mu_off,Nu_off)
               Enddo
            Enddo
         Enddo 

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The symmetry packed AO integrals"
      call output(W_out,1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
#endif 

      Return
      End
