      Subroutine Print_basis_info(Nshells,Ncenters,Nalpha,Npcoef,
     &                            Naobasis,Nbasis,Lnp1,Lnpo,Itfct,
     &                            Ncfps,Npfps,Iqmstr,Nucstr,
     &                            Nrcstr,Jcostr,Jstrt,Jrs,Angmom,
     &                            Ioff_alpha,Ioff_Pcoef,Alpha,Pcoef,
     &                            Shell_coords,Atom_coords,Maxjco,
     &                            Iecp)

      Implicit Integer(A-Z)

      Double Precision Alpha(Nalpha),Pcoef(Npcoef)
      Double Precision Atom_Coords(3,Ncenters)
      Double Precision Shell_Coords(3,NShells)

      Dimension Ncfps(Nshells),Npfps(Nshells),Iqmstr(Ncenters)
      Dimension Jcostr(Ncenters,Maxjco)
      Dimension Nucstr(Ncenters,Maxjco,Maxjco)
      Dimension Nrcstr(Ncenters,Maxjco,Maxjco)
      Dimension Jstrt(Nshells)
      Dimension Angmom(Nshells)
      Dimension Jrs(Nshells)
      Dimension Ioff_pcoef(Nshells)
      Dimension Ioff_alpha(Nshells)
     
      Write(6,*) 
      Write(6,"(a,5(1x,I4))") " Natoms,Nshells,Itfct,Lnp1,Lnpo:",
     &                           Ncenters,Nshells,Itfct,Lnp1,Lnpo
      Write(6,"(a,4(1x,I4))") " Nbasis,Naobasis,Nalpha,Npcoef :",
     &                           Nbasis,Naobasis,Nalpha,NPcoef
      Write(6,*)
      Write(6,"(a)") "Angular momentum"
      Write(6,"(6(1x,I4))") (Angmom(i), i=1, Nshells)
      Write(6,"(a)") "The number of con. function per shell"
      Write(6,"(6(1x,I4))") (nCFpS(i), i=1, nshells)
      Write(6,"(a)") "The number of prim. function per shell"
      Write(6,"(6(1x,I4))") (nPFpS(i), i=1, nshells)
      Write(6,"(a)") "The Ioff_alpha"
      Write(6,"(6(1x,I4))") (Ioff_alpha(i), i=1, Nshells)
      Write(6,"(a)") "The Ioff_pcoef"
      Write(6,"(6(1x,I4))") (Ioff_Pcoef(i), i=1, Nshells)
C 
      If (Iecp .Gt. 0) Then
      Write(6,"(a)") "The Iqmstr"
      Write(6,"(6(1x,I4))") (Iqmstr(i), i=1, ncenters)
      Write(6,"(a)") "The jcostr"
      Do i=1, Ncenters
         Write(6,"(4(1x,I4))") (jcostr(i,j),j=1,Nshells)
      Enddo
      Write(6,"(a)") "The Nucstr"
      Do i=1, Ncenters
         Write(6,"(10(1x,I4))") ((Nucstr(i,j,k),j=1,Nshells),
     &                            k=1,Nshells)
      Enddo
      Write(6,"(a)") "The Nrcstr"
      Do i=1, Ncenters
         Write(6,"(10(1x,I4))") ((Nrcstr(i,j,k),j=1,Nshells),
     &                            k=1,Nshells)
      Enddo
      Write(6,"(a)") "The jstrt"
      Write(6,"(6(1x,I4))") (Jstrt(i), i=1, Nshells)
      Write(6,"(a)") "The jrs"
      Write(6,"(6(1x,I4))") (Jrs(i), i=1, Nshells)
      Write(6,"(a)") "The jcostr"
      Do i=1, Ncenters
         Write(6,"(4(1x,I4))") (jcostr(i,j),j=1,4)
      Enddo
      Endif 

      Write(6,*)
      Write(6,"(a)") "The exponents"
      Write(6,"(6(1x,F12.8))") (Alpha(i), i=1, Nalpha)
      Write(6,*)
      Write(6,"(a)") "The contraction coefs."
      Write(6,"(6(1x,F12.8))") (Pcoef(i), i=1, Npcoef)

      Write(6,*)
      Write(6,*) "Atom Coordinates of each shell"
      Do Ishell =1, Nshells
         Write(6, "(3(1x,F12.8))") (Shell_Coords(I,Ishell), i=1,3)
      Enddo

      Write(6,*)
      Write(6,*) "Atom Coordinates"
      Do Iatom =1, Ncenters 
         Write(6, "(3(1x,F12.8))") (Atom_Coords(I,Iatom), i=1,3)
      Enddo

      Return
      End 

  
