      Subroutine Pccd_check_hess(Work,Maxcor)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Logical Print

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Irrepx = 1
      Print  = .False.
      Write(6,*)
      List_h = 213
      Nsize = Idsymsz(Irrepx,21,21)
      Call Getall(Work,Nsize,Irrepx,List_h)
      Call checksum("H(Ij,Kl):",Work,Nsize)
      If (Print) Call Pccd_print_all(Work,Nsize,"ijkl",Irrepx)
    
      Print  = .False.
      List_h = 234
      Nsize = Idsymsz(Irrepx,19,19)
      Call Getall(Work,Nsize,Irrepx,List_h)
      Call checksum("H(Ab,cD):",Work,Nsize)
      If (Print) Call Pccd_print_all(Work,Nsize,"abcd",Irrepx)

      List_h = 217
      Nsize = Idsymsz(Irrepx,14,15)
      Call Getall(Work,Nsize,Irrepx,List_h)
      Call checksum("H(Ij,Ab):",Work,Nsize)
   
      List_h = 216
      Nsize = Idsymsz(Irrepx,15,14)
      Call Getall(Work,Nsize,Irrepx,List_h)
      Call checksum("H(Ab,iJ):",Work,Nsize)

C These are not needed for orbital optimizations
#ifdef _NOSKIP
      List_h = 210
      Nsize = Idsymsz(Irrepx,14,18)
      Call Getall(Work,Nsize,Irrepx,List_h)
      Call checksum("H(Ij,Ka):",Work,Nsize)

      List_h = 211
      Nsize = Idsymsz(Irrepx,18,14)
      Call Getall(Work,Nsize,Irrepx,List_h)
      Call checksum("H(Ia,Jk):",Work,Nsize)

      List_h = 230
      Nsize = Idsymsz(Irrepx,15,11)
      Call Getall(Work,Nsize,Irrepx,List_h)
      Call checksum("H(Ab,Ci):",Work,Nsize)

      List_h = 231
      Nsize = Idsymsz(Irrepx,11,15)
      Call Getall(Work,Nsize,Irrepx,List_h)
      Call checksum("H(Ci,Ab):",Work,Nsize)

      List_h = 223
      Nsize = Idsymsz(Irrepx,9,9)
      Call Getall(Work,Nsize,Irrepx,List_h)
      Call checksum("H(AI,BJ):",Work,Nsize)

      List_h = 225
      Nsize = Idsymsz(Irrepx,11,11)
      Call Getall(Work,Nsize,Irrepx,List_h)
      Call checksum("H(Ai,Bj):",Work,Nsize)

      List_h = 218
      Nsize = Idsymsz(Irrepx,18,11)
      Call Getall(Work,Nsize,Irrepx,List_h)
      Call checksum("H(Ai,Bj):",Work,Nsize)
#endif 

      Return 
      End 
