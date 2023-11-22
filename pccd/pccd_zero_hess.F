      Subroutine Pccd_zero_hess(Work,Maxcor)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Logical Print

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Irrepx = 1
      Print  = .True.
      Write(6,*)
      List_h = 213
      Nsize = Idsymsz(Irrepx,14,14)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)
    
      Print  = .False.
      List_h = 234
      Nsize = Idsymsz(Irrepx,15,15)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      List_h = 217
      Nsize = Idsymsz(Irrepx,14,15)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)
   
      List_h = 216
      Nsize = Idsymsz(Irrepx,15,14)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      List_h = 210
      Nsize = Idsymsz(Irrepx,14,18)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      List_h = 211
      Nsize = Idsymsz(Irrepx,18,14)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      List_h = 230
      Nsize = Idsymsz(Irrepx,15,11)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      List_h = 231
      Nsize = Idsymsz(Irrepx,11,15)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      Return 
      End 
