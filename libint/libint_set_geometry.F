      Subroutine Libint_set_geometry(Charge,Centers,Iatom,Coord,
     &                               Atom_charge,Ncenters) 

      Implicit None

      integer Iatom
      integer I,Ncenters

      Double precision Coord(3)
      Double precision Atom_charge 
      Double precision Charge(Ncenters)
      Double precision Centers(Ncenters,3)

      Charge(Iatom) = Atom_charge
      do i = 1, 3
        Centers(Iatom,i) = Coord(i)
      enddo

      Return
      End

