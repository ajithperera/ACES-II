      Subroutine Pccd_gossip(Grd,Grd_stata,Grd_statb,Sgrad_stata,
     +                       Sgrad_statb,Nbas,Ispin)

      Implicit Double Precision(A-H,O-Z)
      Logical Symmetry 

      Dimension Grd_stata(6)
      Dimension Grd_statb(6)
      Dimension SGrd_stata(6)
      Dimension SGrd_statb(6)

#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"

      Common /Symm/Symmetry

      If (Ispin .Eq. 1) Then
          Call Pccd_vstat(Grd,Sgrd_stata,Nbas*Nbas)
          Call Putrec(20,"JOBARC","SGRSTATA",6*Iintfp,Sgrd_stata)
      Elseif (Ispin .Eq. 2) Then
          Call Pccd_vstat(Grd,Sgrd_statb,Nbas*Nbas)
          Call Putrec(20,"JOBARC","SGRSTATB",6*Iintfp,Sgrd_statb)
      Endif

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Scalled gradients (K=-H^-1g)"
      Call output(Grd,1,Nbas,1,Nbas,Nbas,Nbas,1)
      Write(6,*)
      If (Ispin .Eq. 1) Write(6,"(2a)") " The statistics of the scaled",
     +                  " alpha orbital rotation gradient matrix"
      If (Ispin .Eq. 2) Write(6,"(2a)") " The statistics of the scaled",
     +                  " beta  orbital rotation gradient matrix"
      Write(6,"(1x,2a)") "--------------------------------------------",
     +                   "-------"
      Write(6,"(5x,a,5xa,5xa)") "Minimum grad.", "Maximum grad.",
     +                           "RMS grad"
      If (Ispin .Eq. 1) Then
          Write(6,*)
          Write(6,"(3(5x,E12.6))") SGrd_stata(3), SGrd_stata(4),
     +                             SGrd_stata(5)
      Elseif (Ispin .Eq. 2) Then
          Write(6,*)
          Write(6,"(3(5x,E12.6))") SGrd_statb(3), SGrd_statb(4),
     +                             SGrd_statb(5)
      Endif
      Write(6,*)
      Write(6,"(1x,2a)") "--------------------------------------------",
     +                   "-------"
#endif

      Return
      End 
