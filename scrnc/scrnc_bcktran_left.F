
      Subroutine Scrnc_bcktran_left(JJ_respn_oo,JJ_respn_xx,
     +                              Eveca,JJ_respn_xo,
     +                              Nij_rows,Nxj_rows,Nxx_rows,
     +                              Nao,Nmo,Irrep_l)

      Implicit Double Precision (A-H,O-Z)

      Double Precision JJ_respn_oo(Nij_rows)
      Double Precision JJ_respn_xo(Nxj_rows)
      Double Precision JJ_respn_xx(Nxx_rows)
      Dimension Eveca(Nao,Nmo)
      Double Precision Null

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)
#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "sym.com"
#include "flags2.com"
#include "flags.com"
#include "parcc.com"
#include "info.com"
    
      Data Null,One /0.0D0,1.0D0/

      Call Getrec(20,"JOBARC","SCFEVECA",Nao*Nmo*Iintfp,Eveca)

      Ioff_ji = 1
      Ioff_xx = 1
 
      Do Irrep_i = 1, Nirrep 
         Irrep_j = Dirprd(Irrep_i,Irrep_l)

         Nsum = Pop(Irrep_j,1)
         Nrow = Iaopop(Irrep_j)
         Ncol = Pop(Irrep_i,1)
        
         Ioff_x = Ioffao(Irrep_j)
         Ioff_j = Pop(Irrep_j,1)

         If (Min(Nrow,Ncol,Nsum) .NE. 0) Then
            Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Eveca(Ioff_x,Ioff_j),
     +                  Nao,JJ_respn_oo(Ioff_ji),Nsum,Null,
     +                  JJ_respn_xo,Nrow)
         Else
            Call Dzero(JJ_respn_xo,Nrow*Ncol)
         Endif 
         
         Nsum = Pop(Irrep_j,1)
         Ncol = Iaopop(irrep_j)
         Nao2 = Ncol

         Ioff_x = Ioffao(Irrep_i)
         Ioff_i = Pop(Irrep_i,1)

         If (Min(Nrow,Ncol,Nsum) .NE. 0) Then
            Call Xgemm("N","T",Nrow,Ncol,Nsum,One,JJ_respn_xo,Nrow,
     +                  Eveca(Ioff_x,Ioff_i),Nao,Null,
     +                  JJ_respn_xx(ioff_xx),Nrow)
         Else
            Call Dzero(JJ_respn_xx(ioff_xx),Nrow*Ncol)
         Endif 

         Ioff_ji = Ioff_ji + Nsum
         Ioff_xx = Ioff_xx + Nrow * Ncol

       Enddo

       Return
       End
       
