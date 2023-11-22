      Subroutine Pccd_resort_phph_2pphh(Work,Maxcor,List_s,List_t)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D00,1/
      Print*, "in",List_s,List_t

      Ispin   = Ione
      Irrepx  = Ione
      Nsize_s = Idsymsz(Irrepx,11,11)
      Nsize_t = Idsymsz(Irrepx,15,14)

      I000 = Ione
      I010 = I000 + Nsize_s
      I020 = I010 + Nsize_t
      Iend = I020 + Nocco(Ispin)*Nvrto(Ispin)
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_resort_phph_2pphh",
     +                                   Iend,Maxcor)
      Call Getall(Work(I000),Nzise_s,Irrepx,List_s)
      Call SStgen(Work(I000),Work(I010),Nsize_s,Vrt(1,Ispin),
     +             Pop(1,Ispin),Vrt(1,Ispin),Pop(1,Ispin),Work(I020),
     +             Irrepx,"1324")
      Call Getall(Work(I010),Nzise_t,Irrepx,List_t)

      Return
      End
       
     
