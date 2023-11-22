#include  "flags.h"

      Program Dkh

      Implicit Double Precision (A-H, O-Z)

      Logical Spherical, Cartesian
      Logical Contract,Symmetry

#include "icore.com"
#include "istart.com"
#include "ecplims.par"
#include "cvwecp.h"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"

      Call aces_init(icore, i0,icrsiz, iuhf, .true.)

      Contract  = (IFLAGS2(h_IFLAGS2_uncontract) .EQ.0)
      Symmetry  = (IFLAGS(h_IOPPAR_symmetry)  .EQ. 2) 
      Spherical = (IFLAGS(h_IFLAGS_spherical) .EQ. 1) 
      If (Spherical) Cartesian  = .False.

      If (Symmetry .and. Spherical .and. Contract) Then
         Write(6,"(20x,a)") "      -----Warning-------"
         Write(6,"(a,a,a)") " When contracted spherical harmonic",
     +                      " basis set with Abelian symmetry"
         Write(6,"(a,a)")   " is requested the DKH transformtion",
     +                      " is carried out using Cartesian basis."
         Write(6,"(a,a)")   " As a result sym=on or off with",
     +                      " the spherical harmonic basis sets",
     +                      " can give slightly different results."
      Endif

C If no-symmetry, spherical or Cartesian much easier to handle
 
      If (Spherical) Then 
          If (symmetry) then
             Call Dkh__sph_sym(Icore(I0),Icrsiz,Iuhf,Symmetry,
     +                         Contract)
          Else
             Call Dkh__sph_nosym(Icore(I0),Icrsiz,Iuhf,Spherical,
     +                           Cartesian,Symmetry,Contract)
          Endif 
      Else
          Call Dkh__cart(Icore(I0),Icrsiz,Iuhf,Spherical,Cartesian,
     +                   Symmetry,Contract)
      Endif 
    
      Call Aces_fin

      Stop
      End
       
      
       

