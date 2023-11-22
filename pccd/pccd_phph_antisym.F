      Subroutine Pccd_phph_antisym(Work,Maxcor,List_h1,List_h2,List_h3)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D0,1/

      Irrepx   = Ione
      Nsize_h1 = Idsymsz(Irrepx,9,9)
      Nsize_h2 = Idsymsz(Irrepx,18,11)
      Nsize_h3 = Idsymsz(Irrepx,11,11)

      Icheck = Nsize_h1+Nsize_h2+Nsize_h3 
      If (Icheck .Ne. 3*Nsize_h1 .Or. Icheck .Ne. 3*Nsize_h2
     +    .Or. Icheck .Ne. 3*Nsize_h3) Then
          Write(6,"(2a)") " Internal inconsistency occured. UHF",
     +                    " reference is not supported!"
          Call Errex 
      Endif 

      I000 = Ione
      I010 = I000 + Nsize_h2
      I020 = I010 + Nsize_h3
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_phph_antisym",
     +                                   Iend,Maxcor)
      
      Call Getall(Work(I000),Nsize_h2,Irrepx,List_h2)
      Call Getall(Work(I010),Nsize_h3,Irrepx,List_h3)
      Call Daxpy(Nsize_h3,Onem,Work(I010),1,Work(I000),1)
      Call Putall(Work(I000),Nsize_h1,Irrepx,List_h1)

      Return
      End 
      
