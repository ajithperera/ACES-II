      Subroutine Vibavg_lvl0_shifts_1der(Prop_derv1,Cfc_contr,Omega,
     +                                   Nreals,Nmodes,Ndata,Icol,Avg)
      Implicit Double Precision(A-H,O-Z)

      Dimension Prop_derv1(Nreals,Nmodes)
      Dimension Cfc_contr(Nmodes,Ndata)
      Dimension Omega(Nmodes)
      Dimension Avg(Nreals,Ndata)
      Logical Jcc,Sigma 

      Data Done/1.0D0/

      Do K = 1, Ndata
        Do J = 1, Nreals
           Do Imode = 1, Nmodes
              Frq = Done/Omega(Imode+Icol)
              Avg(J,K) = Avg(J,K) - Prop_derv1(J,Imode)*
     +                   Cfc_contr(Imode,K)*Frq
            Enddo
        Enddo
      Enddo 

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(2a)") " Contribution from the first-derivative term",
     +                " to NMR shielding constant"
      do i=1,Ndata
      write(6,"(5(1x,F15.8))") (Avg(j,i),j=1,Nreals)
      enddo
#endif

      Return
      End 
     
