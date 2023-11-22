      Subroutine Get_orbsym(Type,Ip1,Ip2,Ip3,Ip4,Imap,Isym,Nbasis,
     +                      Ip1sym,Ip2sym,Ip3sym,Ip4sym,Index1,
     +                      Index2,Index3,Index4,One_pair)

      Implicit Double Precision(A-H,O-Z)
      Character*2 Type
      Logical One_pair
      Dimension Imap(Nbasis)
      Dimension Isym(Nbasis)

#include "info.com"

      One_pair = (Ip3 .eq. 0 .or. Ip4 .eq. 0)
      If (Type .Eq. "OO") Then
         Ip1 = Ip1 + Nocco(1) + 1
         Ip2 = Ip2 + Nocco(1) + 1
         Ip1sym = Isym(Ip1)
         Ip2sym = Isym(Ip2)
         Index1 = Imap(Ip1)
         Index2 = Imap(Ip2)
         If (.not. one_pair) Then
            Ip3 = Ip3 + Nocco(1) + 1
            Ip4 = Ip4 + Nocco(1) + 1
            Ip3sym = Isym(Ip3)
            Ip4sym = Isym(Ip4)
            Index3 = Imap(Ip3)
            Index4 = Imap(Ip4)
         Endif 
      Endif 

      If (Type .Eq. "VV") Then
         Ip1 = Ip1 + Nocco(1)
         Ip2 = Ip2 + Nocco(1)
         Ip1sym = Isym(Ip1)
         Ip2sym = Isym(Ip2)
         Index1 = Imap(Ip1)
         Index2 = Imap(Ip2)
         If (.not. one_pair) Then
            Ip3 = Ip3 + Nocco(1)
            Ip4 = Ip4 + Nocco(1)
            Ip3sym = Isym(Ip3)
            Ip4sym = Isym(Ip4)
            Index3 = Imap(Ip3)
            Index4 = Imap(Ip4)
         Endif 
      Endif 

      If (Type .Eq. "OV") Then
         Ip1 = Ip1 + Nocco(1) + 1
         Ip2 = Ip2 + Nocco(1) 
         Ip1sym = Isym(Ip1)
         Ip2sym = Isym(Ip2)
         Index1 = Imap(Ip1)
         Index2 = Imap(Ip2)
         If (.not. one_pair) Then
            Ip3 = Ip3 + Nocco(1) + 1
            Ip4 = Ip4 + Nocco(1) 
            Ip3sym = Isym(Ip3)
            Ip4sym = Isym(Ip4)
            Index3 = Imap(Ip3)
            Index4 = Imap(Ip4)
         Endif 
      Endif
      If (Type .Eq. "VO") Then
         Ip1 = Ip1 + Nocco(1)
         Ip2 = Ip2 + Nocco(1) + 1
         Ip1sym = Isym(Ip1)
         Ip2sym = Isym(Ip2)
         Index1 = Imap(Ip1)
         Index2 = Imap(Ip2)
         If (.not. one_pair) Then
            Ip3 = Ip3 + Nocco(1)
            Ip4 = Ip4 + Nocco(1) + 1
            Ip3sym = Isym(Ip3)
            Ip4sym = Isym(Ip4)
            Index3 = Imap(Ip3)
            Index4 = Imap(Ip4)
         Endif 
      Endif 

      Return
      End 

  
