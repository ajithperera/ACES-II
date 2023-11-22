#include "flags.h"
#include "disttype.h"

      Subroutine Scrnc_newsize(Work,Memleft,Iuhf,Irrepx,Nsizec)

      Implicit Double Precision (A-H, O-Z)

      Dimension Work(Memleft)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "flags.com"

      Nsizec = Irpdpd(Irrepx,9)
      If (Iuhf .EQ. 0) Then
         Nsizec = Nsizec  + Idsymsz(Irrepx,13,14)
      Else
        Nsizec = Nsizec + Irpdpd(Irrepx,10)
        Nsizec = Nsizec + Idsymsz(Irrepx,1,3)
        Nsizec = Nsizec + Idsymsz(Irrepx,2,4)
        Nsizec = Nsizec + Idsymsz(Irrepx,13,14)
      Endif

      Return
      End
