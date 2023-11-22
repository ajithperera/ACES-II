      Subroutine Form_1dshift_deriv_eq(Plus,Mins,Equi,Diff,Nreals,Nmode,
     +                                 Step_size,Idl,Nfdiffs,Imode,
     +                                 Itype)

      Implicit Double Precision(A-H,O-Z)
      Double precision Mins
      
      Dimension Plus(Nreals,Nfdiffs),Mins(Nreals,Nfdiffs)
      Dimension Equi(Nreals) 
      Dimension Diff(Nreals,Nmode)
      Data Half,One,Two /0.50D0,1.0D0,2.0D0/

      Step = One/(Step_size*Step_size)

      Do I = 1, Nreals 
         D2  = Plus(I,1) + Mins(I,1) -Two*Equi(I)
         Diff(I,Imode) = D2*Step
      Enddo

#ifdef _DEBUG_LVL0
      If (Itype.Eq.1) Write(6,"(a,a)") "The second derivative",
     +                                 " of NMR shifts."
      Write(6,*)
      Write(6,"(6(1x,F15.8))") (Diff(i,imode),i=1,Nreals)
      Write(6,*)
#endif

      Return
      End 
