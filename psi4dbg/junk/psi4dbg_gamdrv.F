      Subroutine Psi4dbg_gamdrv(Work,Maxcor)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Logical Print

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      CALL PSI4DBG_FORMV1(ICORE,MAXCOR,IUHF)
      CALL PSI4DBG_FORMX(ICORE,MAXCOR,IUHF)

#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(a)") "Checksums of G(ab,cd)"
#endif
      CALL Psi4dbg_GAMMA2(ICORE,MAXCOR,IUHF,0)

#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(a)") "Checksums of G(ij,kl)"
#endif
C      CALL Psi4dbg_GAMMA3(ICORE,MAXCOR,IUHF)

      Irrepx = 1
      Write(6,*)
      List_h = 116
      Nsize = Idsymsz(Irrepx,13,14)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      List_h = 113
      Nsize = Idsymsz(Irrepx,14,14)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)
    
      List_h = 133
      Nsize = Idsymsz(Irrepx,13,13)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      List_h = 123
      Nsize = Idsymsz(Irrepx,9,9)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      List_h = 125
      Nsize = Idsymsz(Irrepx,11,11)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      List_h = 118
      Nsize = Idsymsz(Irrepx,9,10)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      Return 
      End 
