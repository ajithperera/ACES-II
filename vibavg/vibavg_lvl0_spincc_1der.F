      Subroutine Vibavg_lvl0_spincc_1der(Prop_derv1,Cfc_contr,Omega,
     +                                   Nreals,Nmodes,Ndata,Icol,Avg)

      Implicit Double Precision(A-H,O-Z)

      Dimension Prop_derv1(Nreals,Nreals,Nmodes)
      Dimension Cfc_contr(Nmodes,Ndata)
      Dimension Omega(Nmodes)
      Dimension Avg(Nreals,Nreals,Ndata)

#include "const.h"
      Data Done,Tensq/1.0D0,100.D0/

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " The scaled property first derivatives"
      do i=1,nmodes
      call output(Prop_derv1(1,1,i),1,Nreals,1,Nreals,Nreals,Nreals,1)
      enddo
#endif 

      Do K = 1, Ndata
         Do Jcc = 1, Nreals 
            Do Icc = Nreals, Jcc+1, -1
               Do Imode = 1, Nmodes 
                  Frq = Done/Omega(Imode+Icol)
                  Avg(Icc,Jcc,K) = Avg(Icc,Jcc,K) -
     +                             Cfc_contr(Imode,K)*Frq*
     +                             Prop_derv1(Icc,Jcc,Imode)
               Enddo
            Enddo
         Enddo 
      Enddo 

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(2a)") " Contribution from the first-derivative term",
     +                " to NMR spin-spin coupling constant"
      do i=1,Ndata
      call output(Avg(1,1,i),1,Nreals,1,Nreals,Nreals,Nreals,1)
      enddo 
#endif 

      Return
      End 
