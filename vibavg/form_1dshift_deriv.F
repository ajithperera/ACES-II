      Subroutine Form_1dshift_deriv(Plus,Mins,Diff,Nreals,Nmode,
     +                             Step_size,Idl,Nfdiffs,Imode,
     +                             Itype)

      Implicit Double Precision(A-H,O-Z)
      Double precision Mins
      
      Dimension Plus(Nreals,Nfdiffs),Mins(Nreals,Nfdiffs)
      Dimension Diff(Nreals,Nmode)
      Data Half /0.50D0/

      Step = Half/Step_size

C First derivatives of the property

      If (Idl .Eq. 1) Then
         Do I = 1, Nreals
            D    = Plus(I,1) - Mins(I,1)
            Diff(I,Imode) = D*Step
         Enddo

C Second derivatives  of the property

      Elseif (Idl .Eq. 2) Then
         Do I = 1, Nreals 
            D1p  = Plus(I,2) - Plus(I,1)
            D1p  = D1p*step
            D1m  = Mins(I,2) - Mins(I,1)
            D1m  = D1m*Step
            D2   = D1p - D1m
            Ioff = Ioff + 1
            Diff(I,Imode) = D2*Step
         Enddo
      Endif

#ifdef _DEBUG_LVL0
      Write(6,*)
      If (Idl .Eq. 1) Then
         If (Itype.Eq.1) Write(6,"(a,a)") "The first derivative",
     +                                  " of NMR shifts."
         write(6,*)
         Write(6,"(5(1x,F15.8))") (Diff(i,imode),i=1,Nreals)
      Elseif (Idl .Eq. 2) Then
         If (Itype.Eq.1) Write(6,"(a,a)") "The second derivative",
     +                                  " of NMR shifts."
         write(6,*)
         Write(6,"(5(1x,F15.8))") (Diff(i,imode),i=1,Nreals)
      Endif 
      Write(6,*)
#endif

      Return
      End 
