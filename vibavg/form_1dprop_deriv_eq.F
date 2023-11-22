      Subroutine Form_1dprop_deriv_eq(Plus,Mins,Equi,Diff,Nreals,Nmode,
     +                                Step_size,Idl,Nfdiffs,Imode,
     +                                Itype)

      Implicit Double Precision(A-H,O-Z)
      Double precision Mins
      
      Dimension Plus(Nreals,3,Nfdiffs),Mins(Nreals,3,Nfdiffs)
      Dimension Equi(Nreals,3) 
      Dimension Diff(3*Nreals,Nmode)
      Data Half,One,Two /0.50D0,1.0D0,2.0D0/

      Step = One/(Step_size*Step_size)

      Do I = 1, Nreals 
         Do J = 1, 3
            D2  = Plus(I,J,1) + Mins(I,J,1) -Two*Equi(I,J)
            Diff(I,Imode) = D2*Step
         Enddo 
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,*)
      If (Itype.Eq.1) Write(6,"(a,a)") "The second derivative",
     +                                 " of NMR shifts."
      If (Itype.Eq.2) Write(6,"(a,a)") "The second derivative",
     +                                 " of gradients."
         Write(6,"(6(1x,F12.6))") (Diff(i,imode),i=1,3*Nreals)
      Write(6,*)
#endif

      Return
      End 
