      Subroutine Vibavg_lvl0_spincc_2der(Prop_derv2,Exponent,Nreals,
     +                                  Nmodes,Ndata,Avg)

      Implicit Double Precision(A-H,O-Z)

      Dimension Prop_derv2(Nreals,Nreals,Nmodes)
      Dimension Exponent(Nmodes,Ndata)
      Dimension Avg(Nreals,Nreals,Ndata)

      Data Done,Factor14,Factor18/1.0D0,0.25D0,0.125D0/

#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,"(a)") " The scaled property second derivatives"
      do i=1,nmodes
      call output(Prop_derv2(1,1,i),1,Nreals,1,Nreals,Nreals,Nreals,1)
      enddo
#endif
      Do K = 1, Ndata 
         Do J = 1, Nreals 
            Do I = Nreals, J-1, -1
               Do Imode = 1, Nmodes 
                  Avg(I,J,K) = Avg(I,J,K) + Prop_derv2(I,J,Imode)*
     +                         Exponent(Imode,K)
               Enddo
            Enddo
         Enddo 
      Enddo 

#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,"(2a)") " Zeroth-order vibrational and thermal", 
     +                " contribution to NMR spin-spin coupling",
      Write(6,"(a)")  "(before scalling)"
      do i=1,Ndata
      call output(Avg(1,1,i),1,Nreals,1,Nreals,Nreals,Nreals,1)
      enddo
#endif

      Call Dscal(Ndata*Nreals*Nreals,Factor14,Avg,1) 

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(2a)") " Zeroth-order vibrational and thermal", 
     +                " contribution to NMR spin-spin coupling"
      do i=1,Ndata
      call output(Avg(1,1,i),1,Nreals,1,Nreals,Nreals,Nreals,1)
      enddo
#endif
      Return
      End 

     
