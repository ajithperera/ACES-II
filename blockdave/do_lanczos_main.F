      Subroutine Do_lanczos_main(Vecs,Work,Memleft,Iuhf,Irrepx,Nsize,
     +                           Maxmem)

      Implicit Double Precision (A-H,O-Z)

#include "lanczos_vars.com"

      Dimension Work(Memleft),Vecs(Nsize,2*Maxvw+6)

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a)") "---Entered Do_lanczos_main ---"
        Write(6,*)
#endif

      Call  Lanczos_main(Vecs,Work,Memleft,Nsize,Irrepx,Iuhf,
     +                   Maxmem)

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a)") "---Exit Do_lanczos_main ---"
        Write(6,*)
#endif

      Return
      End

      

