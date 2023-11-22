      Subroutine Print_spncc_tensor(Tensor,Iatmchrg,Isotopes,Atmlabel,
     +                              Nreals,Natoms) 

      Implicit Double Precision(A-H,O-Z)

      Dimension Isotopes(Nreals)
      Dimension Iatmchrg(Nreals)
      Dimension Tensor(Nreals,Nreals)

      Character*6 Aname_R(45)
      Character*6 Aname_S(55)
      Character*6 AtmLabel(Nreals)

      Data Aname_R/'[1] H ', '[3]HE ', '[7]LI ', '[9]BE ', '[11]B ',
     +             '[13]C ', '[15]N ', '[17]O ', '[19]F ', '[21]NE',
     +             '[23]NA', '[25]MG', '[27]AL', '[29]SI', '[31]P ',
     +             '[33]S ', '[35]CL', '[39]AR', '[39]K ', '[43]CA',
     +             '[45]SC', '[47]TI', '[51]V ', '[53]CR', '[55]MN',
     +             '[57]FE', '[59]CO', '[61]NI', '[63]CU', '[67]ZN',
     +             '[69]GA', '[73]GE', '[75]AS', '[77]SE', '[79]BR',
     +             '[83]KR', '[14]N ', '[37]CL', '[49]TI', '[2] H ',
     +             '[6]LI ', '[10]B ', '[65]CU', '[71]GA', '[81]BR'/

      DATA Aname_S/ '2 [H] ', '[3] H ', '[7]Be ', '[10]BA', '[14]N ',
     +             '[24]Na', '[26]Al', '[32]P ', '[35]S ', '[36]Cl',
     +             '[41]K ', '[43]K ', '[41]Ca', '[45]Ca', '[47]Ca',
     +             '[43]Sc', '[44]Sc', '[45]Ti', '[49]Ti', '[48]V ',
     +             '[49]V ', '[51]Cr', '[52]Mn', '[53]Mn', '[59]Fe',
     +             '[59]Fe', '[55]Co', '[56]Co', '[57]Co', '[60]Co',
     +             '[63]Ni', '[61]Cu', '[65]Cu', '[67]Cu', '[65]Zn',
     +             '[69]Zn', '[71]Zn', '[67]Ga', '[71]Ga', '[72]Ga',
     +             '[69]Ge', '[71]Ge', '[71]As', '[72]As', '[76]As',
     +             '[77]As', '[73]Se', '[75]Se', '[79]Se', '[76]Br',
     +             '[77]Br', '[81]Br', '[82]Br', '[79]Kr', '[81]Kr'/

      Data Izero /0/

      Do Iatm = 1, Nreals
         If (Isotopes(Iatm) .Eq. 0) Then
             Atmlabel(Iatm) = Aname_R(Iatmchrg(Iatm))
         Endif 
      Enddo

      Call Header('The NMR spin-spin coupling tensor (in Hz)',-1,6)
      Call Prntlo_tensor(6,Tensor,Atmlabel,Nreals,Nreals,Nreals,Nreals)
      
      Return 
      End 

