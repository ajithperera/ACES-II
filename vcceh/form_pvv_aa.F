      Subroutine Form_pvv_aa(T2cc,T2mp,T2,Iwork,Imaxcor,Nab,
     +                       Irrepab,Ipert,Ispin)

      Implicit Double Precision(A-H,O-Z)
      Integer Aend,Bend

#include "syminf.com"
#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "symoff.com"
#include "sympop.com"
#include "active_space.com"

      Dimension Iwork(Imaxcor)
      Dimension T2cc(Nab) 
      Dimension T2mp(Nab) 
      Dimension T2(Nab) 
   
      Data Ione /1/

      Ioff = Ione

      Do Irrep_r = 1, Nirrep

         Irrep_b = Irrep_r
         Irrep_a = Dirprd(Irrep_b,Irrepab)

         Ncol_b  = Vrt(Irrep_b,Ispin)
         Nrow_a  = Vrt(Irrep_a,Ispin)

         Aend    = Pactive_vv(Ipert,Irrep_a,Ispin)
         Bend    = Pactive_vv(Ipert,Irrep_b,Ispin)

         Call Form_vv(T2cc(Ioff),T2mp(Ioff),T2(Ioff),Iwork,Imaxcor,
     +                Nrow_a,Ncol_b,Aend,Bend)

         Ioff = Ioff + Nrow_a*Ncol_b

      Enddo
     
      Return 
      End 
