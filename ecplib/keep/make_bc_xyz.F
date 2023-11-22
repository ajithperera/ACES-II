      Subroutine Make_BC_xyz(Nxyz, CSxyz, TBCxyz)
C
c This routine compute sums; Sum_a{a=0,n_i}{fac(ni,a)C_i^^(n_i-a)}
C where i=x,y,z
C
#include "ecp.par"

      Implicit Double Precision (A-H, O-Z)
    
      Dimension TBCxyz(0:Maxang), CSxyz(0:Maxang)

      II = 0
      TBCxyz(0) = Bicoefs(Nxyz,0)*CSxyz(0)
      If (Nxyz .EQ. 0) Return

      Do I = 1, Nxyz

            Pre_Fac = Bicoefs(Nxyz, I)
            II = II + 1

            TBCxyz(II) = TBCxyz(II-1) + Pre_Fac*CSxyz(II)

      Enddo
    
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "um from Sum_xyz_1"
      Write(6,"(F10.6)") (TBCxyz(I), I=0, Nxyz)
#endif
     
      Return
      End
