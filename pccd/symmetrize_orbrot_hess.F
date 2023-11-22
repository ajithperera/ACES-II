      Subroutine Symmetrize_orbrot_hess(W,WT,WA,Nrow,Ncol,Maxd,Irrepl,
     +                                  Irrepr,Type)
  
      Implicit Double Precision(A-H,O-Z)
      Dimension W(Nrow,Ncol),WT(Ncol,Nrow)
      Dimension Wa(Maxd)
      Character*2 Type

      Data Ione /1/

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Ispin = Ione 
      Do Icol = 1, Ncol
         Ioff  = Ione
         Do Irrep_p = 1, Nirrep 
            Irrep_q = Dirprd(Irrep_p,Irrepl)
            If (Type .Eq. "PP") Then
               Np = Vrt(Irrep_p,Ispin)
               Nq = Vrt(Irrep_q,Ispin)
            Else 
               Np = Pop(Irrep_p,Ispin)
               Nq = Pop(Irrep_q,Ispin)
            Endif
            Call pccd_asymm(W(Ioff,Icol),WA,Np,Nq)
            Call Dcopy(Np*Nq,Wa,1,W(Ioff,Icol),1)
            Ioff = Ioff + Np*Nq
         Enddo  
      Enddo 

      Call Transp(W,WT,Ncol,Nrow)

      Do Irow = 1, Nrow
         Ioff  = Ione
         Do Irrep_p = 1, Nirrep 
            Irrep_q = Dirprd(Irrep_p,Irrepr)
            If (Type .Eq. "PP") Then
               Np = Vrt(Irrep_p,Ispin)
               Nq = Vrt(Irrep_q,Ispin)
            Else 
               Np = Pop(Irrep_p,Ispin)
               Nq = Pop(Irrep_q,Ispin)
            Endif
            Call pccd_asymm(WT(Ioff,Irow),Wa,Np,Nq)
            Call Dcopy(Np*Nq,Wa,1,WT(Ioff,Irow),1)
            Ioff = Ioff + Np*Nq
         Enddo  
      Enddo 

      Call Transp(WT,W,Nrow,Ncol)

      Return
      End
