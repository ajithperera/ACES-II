      Subroutine Tdee_pes_ode(T,Data_in_out,Mubar_dot)

      Implicit Double Precision (A-H,O-Z)
      Double Precision Mubar_dot, Mone
      Logical Real

#include "tdee_vars.com"
#include "icore.com"
#include "istart.com"
#include "machsp.com"

      Dimension Data_in_out(2*Nsize_dummy),Mubar_dot(2*Nsize_dummy)
      Data Mone /-1.0D0/

      Irrepx  = Irrepx_dummy
      Iside   = Iside_dummy
      Iuhf    = Iuhf_dummy 
      Nsize   = Nsize_dummy
      Memleft = Memleft_dummy 

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "Entered - tdee_pes_ode"
      Write(6,"(a,3i2,3i10)") "The global vars: ", Irrepx,Iside,
     +                        Iuhf,Nsize,Memleft
CSSS      call prvecr(Data_in_out,2*Nsize)
CSSS      call checksum("@-tdee_pes_ode;Mubar:",Data_in_out,2*Nsize,s)
#endif

C First do the -i Hbar * Mubar_R and put as the imaginary component 
C of the derivative (Real=.True. means we are working on the real 
C component of Mubar(t)).

      Real = .True. 

      Call Tdee_get_derivatives(Irrepx,Iside,Iuhf,Memleft,Nsize,
     +                          Real,Data_in_out,Mubar_dot(Nsize+1))

C Now do the -i*i Hbar * Mbar_I and put as the real component 
C of the derivative.

      Real = .False.

      Call Tdee_get_derivatives(Irrepx,Iside,Iuhf,Memleft,Nsize,
     +                          Real,Data_in_out(Nsize+1),Mubar_dot)

#ifdef _DEBUG_LVL0
CSSS      call checksum("@-tdee_pes_ode;Mubar_dot:",Mubar_dot,2*Nsize,s)
#endif 

      Return 
      End

 
