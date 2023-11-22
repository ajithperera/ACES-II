      Subroutine checkvv(H,Nbas,Nvrt)

      Implicit Double Precision (A-H,O-Z)
      Dimension H(Nbas,Nbas)
      Dimension Ioffo(8)
      Dimension Ioffv(8)
      Integer An,Bn
   
#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"

      Ioffo(1) = 0
      Ioffv(1) = Nocco(1)

      Do Irrep =2, Nirrep
         Ioffo(Irrep)=Ioffo(Irrep-1)+Pop(irrep-1,1)
         Ioffv(Irrep)=Ioffv(Irrep-1)+Vrt(irrep-1,1)
      Enddo

      Ispin = 1
      e = 0.0D0
      Do Irrep_b = 1, Nirrep
         Irrep_a = Dirprd(Irrep_b,1)
         Do b = 1, Vrt(Irrep_b,Ispin)
            Do a = 1, Vrt(Irrep_a,Ispin)
               An = A + Ioffv(Irrep_a)
               Bn = B + Ioffv(Irrep_b)
               e = e +  (H(an,bn)*H(an,bn))
            Enddo
          Enddo
      Enddo
      Print*, E
      Return
      End 
