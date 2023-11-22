      Subroutine Get_lnmodes(Work, Maxcor, Nrel_atoms, Ntot_atoms)

      Implicit Double Precision (A-H, O-Z)

      Dimension Distance(Nrel_atoms, Nrel_atoms)
      Double Precision Nmodes(3*Nrel_atoms, 3*Nrel_atoms), 
     &                 Coord(3,Nrel_atoms), 
     &                 Tmp_coord(3,Ntot_atoms),
     &                 Center(Nrel_atoms, 3*Nrel_atoms)
      Integer Iatmchrg(Ntot_atoms)

#include "machsp.com"      
 
      NRX = 3*Nrel_atoms
      Call Getrec(20, "JOBARC", "SYM_NMDS", NRX*NRX*IINTFP, Nmodes)
      Call Getrec(20, "JOBARC", 'LINEAR  ', 1, ILinear)
      Call Getrec(20, "JOBARC", "COORD  ", 3*Ntot_atoms,Tmp_Coord)
      Call Getrec(20, 'JOBARC', 'ATOMCHRG', Ntot_atoms, Iatmchrg)

      If (Ilinear .EQ. 1) Then
          Nvibs = Nrx - 5
          NRTs  = 5
      Else
          Nvibs = Nrx - 6
          NRTs  = 6
      Endif 
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The Normal modes"
      Do Imode = 1, Nvibs
         Joff = NRTs*3*Nrel_atoms
         Write(6,"(3(F20.10))") (Nmodes(Joff+Jcord, Imode), 
     &                                  Jcord=1, NRX)
      Enddo
      Write(6,*)
      Write(6,"(a)") "The Cartesian Coordiantes of atoms"
      Do iatom = 1, Ntot_atoms
         Write(6,"(3(F20.10))") (Tmp_Coord(Ixyz, Iatom), Ixyz=1, 3)
      Enddo
      Write(6,*)
      Write(6,"(6(i4))") (Iatmchrg(I), I=1, Ntot_atoms)
#endif
      Jatm = 0
      Do Irel_atom = 1, Ntot_atoms
         If (Iatmchrg(Irel_atom) .NE. 0) Then
             Jatm = Jatm + 1
             Call Dcopy(3, Tmp_Coord(1, IRel_atom), 1, 
     &                  Coord(1, Jatm), 1)
         Endif 
      Enddo 
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The Cartesian Coordiantes of real atoms"
      Do iatom = 1, Nrel_atoms 
         Write(6,"(3(F20.10))") (Coord(Ixyz, Iatom), Ixyz=1, 3)
      Enddo
#endif
C
      Do Imode = 1, Nvibs 
         Joff =  NRTs*3*Nrel_atoms 
         Do Iatm = 1, Nrel_atoms
             Write(6,*) Joff, Imode
             Center(Iatm, Imode) = 
     &                     Nmodes(Joff+1, imode)*Nmodes(Joff+1, imode)
     &                   + Nmodes(Joff+2, imode)*Nmodes(Joff+2, imode)
     &                   + Nmodes(Joff+3, imode)*Nmodes(Joff+3, imode)
         Joff = Joff + 3
         Enddo
      Enddo
C 
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The centers of normal modes"
      Do imode = 1, Nvibs
         Do iatm = 1, Nrel_atoms
            Write(6,"((F20.10))") Center(Iatm, Imode)
         Enddo
      Enddo
#endif

C Build a distance matirix among normal modes.
C
C      DO Iatm = 1, Natoms 
C         Do Jatm 1, Natoms
C            Do Ixyz = 1, 3
C               Do Imodes = 1, 3*Natoms
C                  Distance(Iatm, Jatm, Imodes) = Nmodes(Ia
C      
      Return
      End
