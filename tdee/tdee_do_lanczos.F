      Subroutine Tdee_do_lanczos(Vecs,Work,Memleft,Iuhf,Irrepx,Nsize,
     +                           Maxmem)

      Implicit Double Precision (A-H,O-Z)

#include "lanczos_vars.com"

      Dimension Work(Memleft),Vecs(Nsize,2*Maxvw+6)

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a)") "---Entered Tdee_do_lanczos ---"
        Write(6,*)
#endif

      Call  Tdee_Lanczos_main(Vecs,Work,Memleft,Nsize,Irrepx,Iuhf,
     +                        Maxmem)

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a)") "---Exit Tdee_do_lanczos ---"
        Write(6,*)
#endif

      Return
      End

      

