      Subroutine Form_1dprop_deriv(Plus,Mins,Diff,Nreals,Nmode,
     +                             Step_size,Idl,Nfdiffs,Imode,
     +                             Itype)

      Implicit Double Precision(A-H,O-Z)
      Double precision Mins
      
      Dimension Plus(Nreals,3,Nfdiffs),Mins(Nreals,3,Nfdiffs)
      Dimension Diff(3*Nreals,Nmode)
      Data Half /0.50D0/

      Step = Half/Step_size

C First derivatives of the property

      If (Idl .Eq. 1) Then
         Ioff = 0
         Do I = 1, Nreals
            Do J = 1, 3
               D    = Plus(I,J,1) - Mins(I,J,1)
               Ioff = Ioff + 1
               Diff(Ioff,Imode) = D*Step
            Enddo 
         Enddo

C Second derivatives  of the property

      Elseif (Idl .Eq. 2) Then
         Ioff = 0
         Do I = 1, Nreals 
            Do J = 1, 3
               D1p  = Plus(I,J,2) - Plus(I,J,1)
               D1p  = D1p*step
               D1m  = Mins(I,J,2) - Mins(I,J,1)
               D1m  = D1m*Step
               D2   = D1p - D1m
               Ioff = Ioff + 1
               Diff(Ioff,Imode) = D2*Step
            Enddo 
         Enddo
      Endif

#ifdef _DEBUG_LVL0
      Write(6,*)
      If (Idl .Eq. 1) Then
         If (Itype.Eq.1) Write(6,"(a,a)") "The first derivative",
     +                                  " of NMR shifts."
         If (Itype.Eq.2) Write(6,"(a,a)") "The first derivative",
     +                                  " of gradients."
         write(6,*)
         Write(6,"(5(1x,F15.10))") (Diff(i,imode),i=1,3*Nreals)
      Elseif (Idl .Eq. 2) Then
         If (Itype.Eq.1) Write(6,"(a,a)") "The second derivative",
     +                                  " of NMR shifts."
         If (Itype.Eq.2) Write(6,"(a,a)") "The second derivative",
     +                                  " of gradients."
         write(6,*)
         Write(6,"(5(1x,F15.10))") (Diff(i,imode),i=1,3*Nreals)
      Endif 
      Write(6,*)
#endif

      Return
      End 
