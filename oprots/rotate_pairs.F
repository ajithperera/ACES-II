      Subroutine Rotate_pairs(Scf_mos,Eps,Index1,Index2,Ip1Sym,Ip2sym,
     +                         Nbasis,Type)

      Implicit Double Precision(A-H,O-Z)
      Character*2 Type 

      Dimension Scf_mos(Nbasis,Nbasis)

      Cos_eps = Cos(Eps)
      Sin_eps = Sin(EPs)

C For p' = cos (eps) p + Sin (eps) q (p and q are original MOs. Here the
C the number of rows goes over rhe number of function but this can
C easily change go over only the number of orbitals per irrep. Just a
C lot of work to gain nothing!

      If (Type .EQ. "OO" .or. Type .EQ. "VV" .OR. Type .EQ. "VO") Then
      Do Irow = 1, Nbasis

         Scf_mos(Irow,Index1) = Cos_eps*Scf_mos(Irow,INdex1) +
     +                          Sin_eps*Scf_mos(Irow,INdex2) 

         Scf_mos(Irow,Index2) = Sin_eps*Scf_mos(Irow,INdex1) -
     +                          Cos_eps*Scf_mos(Irow,INdex2) 

      Enddo

      ELse if (Type .EQ. "OV") Then

      Do Irow = 1, Nbasis

         Scf_mos(Irow,Index2) = Cos_eps*Scf_mos(Irow,INdex1) +
     +                          Sin_eps*Scf_mos(Irow,INdex2) 

         Scf_mos(Irow,Index1) =  Sin_eps*Scf_mos(Irow,INdex1) -
     +                           Cos_eps*Scf_mos(Irow,INdex2) 

      Enddo
      Endif 

      Return
      End

