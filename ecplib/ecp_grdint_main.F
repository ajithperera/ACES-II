#include "flags.h"

      Subroutine ECP_grdint_main(Ntotatoms, Dens_fao, Naobfns, Grad_xyz)
C      
      Implicit Double Precision(A-H, O-Z)

#include "ecp.par"
#include "ecp.com"
#include "flags.com"
#include "flags2.com"

      Dimension Cint(Maxmem), Dens_fao(Naobfns,Naobfns),
     &          Grad_xyz(3, Ntotatoms), Coord(3,Mxatms), IGenby(Mxatms)

      Call Dzero(Cint, Maxmem)
      Call Getrec(20, 'JOBARC', 'COORD   ', 3*Ntotatoms, Coord)
      Call Symeqv(Ntotatoms, IGenby)
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The symmetry equivalent atoms"
      Write(6,"(6(1x,I3))") (IGenby(i), i=1, Ntotatoms)
      Write(6,*)
      Do j=1, Ntotatoms
      Write(6, "(3(1x,F10.6))") (Coord(i,j), i=1, 3)
      Enddo
      Write(6,*)
#endif

      Call ECP_grdint_driver(Ntotatoms, IGenby, Coord, Cint,
     &                       Dens_fao, Naobfns, Grad_xyz, .True.)
      
      Return
      End
