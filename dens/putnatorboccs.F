      Subroutine Putnatorboccs(Occ_nums,Nbas,Type)
   
      Implicit Double Precision (A-H,O-Z)

      Double Precision Occ_nums(Nbas,Nbas) 
      Character*5 Type
      Character*1 Blank

      Logical Natorbs_occ_exist

      Blank =""
      Inquire(File="NATORB_OCCS",EXIST=Natorbs_occ_exist)
      If (.NOT. Natorbs_occ_exist) Then
         Open(Unit=99,File="NATORB_OCCS",Status="New",Access=
     +       "Sequential",Form="Formatted")
         If (Type .EQ. "Alpha") Then
            Do i = 1, Nbas
               Write(99,"(F7.4)")  Occ_nums(I,I)
            Enddo
         Else
            Write(6,"(2a)") " NATORB_OCCS does not exist. First call",
     +                      " this requires that the third argument" 
            Write(6,"(a)") " is Alpha."
            Call Errex
         Endif 
         
      Else if (Type .EQ. "Beta ") Then

         Write(99,"(a)") BLANK
         Do i = 1, Nbas
            Write(99,"(F7.4)")  Occ_nums(I,I)
         Enddo
      Endif 

      Return
      End 
