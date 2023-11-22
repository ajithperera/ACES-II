      Subroutine Pccd_setg(Work,Maxcor)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Irrepx = 1

#ifdef _NOSKIP
      List_h = 116
      Nsize = Idsymsz(Irrepx,15,14)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)
#endif 

      List_h = 123
      Nsize = Idsymsz(Irrepx,9,9)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      List_h = 125
      Nsize = Idsymsz(Irrepx,11,11)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      List_h = 118
      Nsize = Idsymsz(Irrepx,18,11)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)
      List_h = 113
      Nsize = Idsymsz(Irrepx,14,14)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      List_h = 133
      Nsize = Idsymsz(Irrepx,15,15)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      Return 
      End 
