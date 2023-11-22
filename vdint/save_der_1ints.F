      Subroutine Save_der_1ints(Maxlop,Kfdao,Work,Lwork,Save,Retrieve,
     +                          Kinetic,Attract)

      Implicit Double Precision (A-H, O-Z)
      Dimension  Work(Lwork),Kfdao(8)
      Logical Save,Retrieve,Kinetic,Attract

      Common /Pert/Ntpert,Npert(8),Ipert(8),Ixpert,Iypert,Izpert,
     +             Iyzpert,Ixzpert,Ixypert,Itransx,Itransy,Itransz,
     +             Nucind1
      Common /Lsym/Nlenq(8),Nlent(8)

      If (Save .Or. Retrieve) Then
          Open(unit=11,file="Onehd",form="unformatted",Status="unknown",
     +     access="sequential") 
          Open(unit=12,file="Attrd",form="unformatted",Status="unknown",
     +     access="sequential") 
      Endif 

      If (Save .And. Kinetic) Then
         Do I = 1, Maxlop+1
            Write(Unit=11) Nlent(i)
            Write(Unit=11) (Work(Kfdao(i)+j-1),j=1,Nlent(i))
         Enddo
      Endif

      If (Save .And. Attract) Then
         Do I = 1, Maxlop+1
            Write(Unit=12) Nlent(i)
            Write(Unit=12) (Work(Kfdao(i)+j-1),j=1,Nlent(i))
         Enddo
      Endif

      If (Retrieve .And. Kinetic) Then
         Do I = 1, Maxlop+1
            Read(Unit=11) Nlent(i)
            Read(Unit=11) (Work(Kfdao(i)+j-1),j=1,Nlent(i))
         Enddo
      Do I = 1, Maxlop+1
         Read(Unit=11) Nlent(i)
         Write(6,"(5(1x,F15.10))") (Work(Kfdao(i)+j-1),j=1,Nlent(i))
      Enddo
      Endif
      Return 

      Do Irrep = 1, Maxlop+1
         If (Nlent(Irrep) .Gt. 0) Then
            Mpert=Npert(Irrep)
            Jtransy = Itransy

            If (Ixpert .Eq. Irrep) Then
               Mpert = Mpert-1
               Call Reduce(Work(Kfdao(Irrep)),Nlent(Irrep),Itransx,
     +                     Mpert+1)
               If (Ixpert .EQ. Iypert .and. Itransy .Gt. Itransx)
     +             Jtransy=Jtransy-1
               If (Ixpert .EQ. Izpert .and. Itransz .Gt. Itransx)
     +             Jtransz=Jtransz-1
            Endif 

           If (Iypert .Eq. Irrep) Then
               Mpert = Mpert-1
               Call Reduce(Work(Kfdao(Irrep)),Nlent(Irrep),Jtransy,
     +                     Mpert+1)
               If (Iypert .EQ. Izpert .and. Itransz .Gt. Itransy)
     +             Jtransz=Jtransz-1
            Endif

           If (Izpert .Eq. Irrep) Then
               Mpert = Mpert-1
               Call Reduce(Work(Kfdao(Irrep)),Nlent(Irrep),Jtransz,
     +                     Mpert+1)
            Endif
         Endif  
      Enddo 
    
      Return
      End
