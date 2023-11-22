      Subroutine Get_isodata(Atmmass,Atmlabel,Iatmchrg,Isotopes,
     +                       Name,Iunitis,Nreals)

      Implicit Double Precision(A-H,O-Z)

      Character*5 AtmLabel(Nreals)
      Dimension Atmmass(Nreals)
      Dimension Iatmchrg(Nreals)
      Dimension Isotopes(Nreals)
      Dimension I(Nreals)
      Dimension D(Nreals)
      Character*3 C(Nreals)
      Character*6 Name(Nreals)
      Character*6 Aname_R(45)
      Character*6 Aname_s(55)
      Character*1 V1
      Character*2 V2
      Logical Found 

      Data Half,Inull /0.50d0,0/
      DATA ANAME_R/ '[1] H ', '[3]He ', '[7]Li ', '[9]Be ', '[11]B ',
     &              '[13]C ', '[15]N ', '[17]O ', '[19]F ', '[21]Ne',
     &              '[23]Na', '[25]Mg', '[27]Al', '[29]Si', '[31]P ',
     &              '[33]S ', '[35]Cl', '[39]Ar', '[39]K ', '[43]Ca',
     &              '[45]Sc', '[47]Ti', '[51]V ', '[53]Cr', '[55]Mn',
     &              '[57]Fe', '[59]Co', '[61]Ni', '[63]Cu', '[67]Zn',
     &              '[69]Ga', '[73]Ge', '[75]As', '[77]Se', '[79]Br',
     &              '[83]Kr', '[14]N ', '[37]Cl', '[49]Ti', '[2] H ',
     &              '[6]Li ', '[10]B ', '[65]Cu', '[71]Ga', '[81]Br'/

      DATA ANAME_S/ '2 [H] ', '[3] H ', '[7]Be ', '[10]BA', '[14]N ',
     &              '[24]Na', '[26]Al', '[32]P ', '[35]S ', '[36]Cl',
     &              '[41]K ', '[43]K ', '[41]Ca', '[45]Ca', '[47]Ca',
     &              '[43]Sc', '[44]Sc', '[45]Ti', '[49]Ti', '[48]V ',
     &              '[49]V ', '[51]Cr', '[52]Mn', '[53]Mn', '[59]Fe',
     &              '[59]Fe', '[55]Co', '[56]Co', '[57]Co', '[60]Co',
     &              '[63]Ni', '[61]Cu', '[65]Cu', '[67]Cu', '[65]Zn',
     &              '[69]Zn', '[71]Zn', '[67]Ga', '[71]Ga', '[72]Ga',
     &              '[69]Ge', '[71]Ge', '[71]As', '[72]As', '[76]As',
     &              '[77]As', '[73]Se', '[75]Se', '[79]Se', '[76]Br',
     &              '[77]Br', '[81]Br', '[82]Br', '[79]Kr', '[81]Kr'/


      Write(6,"(15x,a)") "-------------Alert---------------------"
      Write(6,"(5x,2a)") "Atom label, charge and mass is read from the",
     +                   " user given iSOTOPES file"
      Write(6,"(15x,a)") "---------------------------------------"

#ifdef _DEBUG_LVL0
      write(6,"(a)") " At entry Get_isodata"
      Write(6,"(6(1x,I4))") (Iatmchrg(j),j=1,Nreals)
      Write(6,"(6(1x,F10.6))") (Atmmass(j),j=1,Nreals)
      Write(6,"(6(1x,6a))") (Atmlabel(j),j=1,Nreals)
#endif 

      Do Iatm = 1, Nreals
         Read(Iunitis,*)  C(Iatm), I(Iatm), D(Iatm)
         If (I(Iatm)-Iatmchrg(Iatm) .Ne. Inull) Then
             Write(6,10) " Isotopic mass", D(Iatm)," and charge",
     +                     I(Iatm), " is used for atom", Iatm, 
     +                   " and label ", C(iatm)
             Isotopes(Iatm) = 1
         Endif
         Atmmass(Iatm)  = D(Iatm)
         Iatmchrg(Iatm) = I(Iatm)
         AtmLabel(Iatm) = C(Iatm)
      Enddo

   10 Format(1x,a,F12.6,a,1x,I3,a,1x,I2,2a)

      Do Iatm = 1, Nreals
         If (Iatmchrg(Iatm) .Lt. 10) Then
             Write(V1,"(I1)") Iatmchrg(Iatm)
             Name(Iatm) = "["//V1//"] "//AtmLabel(Iatm)
         Elseif  (Iatmchrg(Iatm) .Lt. 100) Then
             Write(V2,"(I2)") Iatmchrg(Iatm)
             Name(Iatm) = "["//V2//"]"//AtmLabel(Iatm)
         Endif 
      Enddo

      Do Iatm = 1, Nreals
         Atmlabel(Iatm) = Name(Iatm)
      Enddo 

#ifdef _DEBUG_LVLM
      write(6,*)
      Write(6,"(a)") " The isotopes array;unmodified"
      Write(6,"(6(1x,I4))") (Isotopes(j),j=1,Nreals)
#endif

      Do Iatm = 1, Nreals
         Icount = 0
         Found  = .False.
         Do while (Found)
            Icount = Icount + 1
            If (Aname_s(Icount) .Eq. Name(Iatm)) Then
               Found = .True.
               Isotopes(Iatm) = Icount
            Elseif (Icount .Eq. 38) Then
               Found = .True.
            Endif
         Enddo
      Enddo 

#ifdef _DEBUG_LVL0
      write(6,*)
      Write(6,"(6(a,1x,6a))") " Istope labels:",(Name(j),j=1,Nreals)
      Write(6,"(a)") " The isotopes array:modified"
      Write(6,"(6(1x,I4))") (Isotopes(j),j=1,Nreals)
      Write(6,"(6(1x,I4))") (Iatmchrg(j),j=1,Nreals)
#endif

      Return
      End


