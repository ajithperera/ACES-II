      Subroutine Form_2dprop_deriv_eq(Plus,Mins,Equi,Diff,Ndim,Nmodes,
     +                                Step_size,Idl,Nfdiffs,Imode,
     +                                Itype)

      Implicit Double Precision(A-H,O-Z)
      Double Precision Mins
      
      Dimension Plus(Ndim,Ndim,Nfdiffs),Mins(Ndim,Ndim,Nfdiffs)
      Dimension Equi(Ndim,Ndim)
      Dimension Diff(Ndim,Ndim,Nmodes)
      Data Half,One,Two /0.50D0,1.0D0,2.0D0/

      Step = One/(Step_size*Step_size)

#ifdef _DEBUG_LVL0
      Write(6,"(2a)") "In Form_2dprop_deriv_eq property matrices."
      Write(6,"(2a)") "Calculating 2nd derivative using equilibrium"
     +                " and plu/minus dispalced values"
      Write(6,"(a)") "Plus"
      call output(plus(1,1,1),1,Ndim,1,Ndim,Ndim,Ndim,1)
      Write(6,"(a)") "Minus"
      call output(mins(1,1,1),1,Ndim,1,Ndim,Ndim,Ndim,1)
      Write(6,"(a)") "Equi"
      call output(Equi(1,1),1,Ndim,1,Ndim,Ndim,Ndim,1)
#endif

C Second derivatives of the property

      Do I = 1, Ndim
         Do J = 1, Ndim
            D2 = Plus(I,J,1) + Mins(I,J,1) - Two*Equi(I,J)
            Diff(I,J,Imode) = D2*Step
         Enddo 
      Enddo 

#ifdef _DEBUG_LVL0
      If (Itype.Eq.1)
     + write(6,"(a,1x,I3)") " The second-derivative of JCC for mode",
     +                        Imode
      call output(Diff(1,1,Imode),1,Ndim,1,Ndim,Ndim,Ndim,1)
#endif 
     
      Return
      End 

      

