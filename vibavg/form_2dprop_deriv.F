      Subroutine Form_2dprop_deriv(Plus,Mins,Diff,Ndim,Nmodes,
     +                             Step_size,Idl,Nfdiffs,Imode,
     +                             Itype)

      Implicit Double Precision(A-H,O-Z)
      Double Precision Mins
      
      Dimension Plus(Ndim,Ndim,Nfdiffs),Mins(Ndim,Ndim,Nfdiffs)
      Dimension Diff(Ndim,Ndim,Nmodes)
      Data Half /0.50D0/

      Step = Half/step_size

#ifdef _DEBUG_LVL0
      Write(6,"(2a)") "In Form_2dprop_deriv property matrices.",
     +                " 2 matrice for first and 4 matrices for"
      Write(6,"(a)")  " second derivatives"
      Write(6,"(a)") " + c+deltax*q"
      call output(plus(1,1,1),1,Ndim,1,Ndim,Ndim,Ndim,1)
      call output(mins(1,1,1),1,Ndim,1,Ndim,Ndim,Ndim,1)
      If (Idl .eq. 2) Then
      Write(6,*) 
      Write(6,"(a)") " + c-deltax*q"
      call output(plus(1,1,2),1,Ndim,1,Ndim,Ndim,Ndim,1)
      call output(mins(1,1,2),1,Ndim,1,Ndim,Ndim,Ndim,1)
      endif 
#endif
C First derivatives of the property 

      If (Idl .Eq. 1) Then
         Do I = 1, Ndim
            Do J = 1, Ndim
               D = Plus(J,I,Idl) - Mins(J,I,Idl)
               Diff(J,I,Imode) = D*Step
            Enddo
         Enddo 

#ifdef _DEBUG_LVL0
      If (Itype .Eq. 1) 
     + write(6,"(a,1x,I3)") " The first-derivative of JCC for mode",
     +                        Imode
      If (Itype .Eq. 2) 
     + write(6,"(2a,1x,I3)") " The first-derivative of Hessian for",
     +                        "  mode", Imode
      call output(Diff(1,1,Imode),1,Ndim,1,Ndim,Ndim,Ndim,1)
#endif 

      Elseif (Idl .Eq. 2) Then

C Second derivatives  of the property

         Do I = 1, Ndim
            Do J = 1, Ndim
               D1p = Plus(I,J,2) - Plus(I,J,1)
               D1p = D1p*step
               D1m = Mins(I,J,2) - Mins(I,J,1)
               D1m = D1m*Step
               D = D1p - D1m
               Diff(I,J,Imode) = D*Step
            Enddo 
         Enddo 
#ifdef _DEBUG_LVL0
      If (Itype.Eq.1)
     + write(6,"(a,1x,I3)") " The second-derivative of JCC for mode",
     +                        Imode
      call output(Diff(1,1,Imode),1,Ndim,1,Ndim,Ndim,Ndim,1)
      If (Itype .Eq. 2) 
     + write(6,"(2a,1x,I3)") " The second-derivative of Hessian for",
     +                        "  mode", Imode
      call output(Diff(1,1,Imode),1,Ndim,1,Ndim,Ndim,Ndim,1)
#endif 
      Endif 
     
      Return
      End 

      

