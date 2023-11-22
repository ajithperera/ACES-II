#include "flags.h"
      Subroutine Get_mbpt2_t2s(T2amp,Work,Maxcor,T2ln,T2ln_aa,T2ln_bb,
     +                         T2ln_ab,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Integer T2ln,T2ln_aa,T2ln_bb,T2ln_ab
      Dimension Work(Maxcor),T2amp(T2ln)

#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "flags.com"
#include "flags2.com"
#include "sympop.com"
#include "info.com"

      Do Ispin = 1, (Iuhf+1)
         Listw = 13 + Ispin
         Listd = 47 + Ispin
         Ibegin = (Ispin-1)*T2ln_aa
         Ioff   =  1
         Do irrep = 1, Nirrep
            Nrow = Irpdpd(Irrep,Ispin)
            Ncol = Irpdpd(Irrep,2+Ispin)
            I000 = 1 
            I010 = I000 + Nrow*Ncol
            Iend = I010 + Nrow*Ncol

            If (Iend .Ge. Maxcor) Call Insmem("get_mbpt_nos",
     +                                         Iend,Maxcor)
            Call Getlst(Work(I000),1,Ncol,2,Irrep,Listw)
            Call Getlst(Work(I010),1,Ncol,2,Irrep,ListD)
            Call vecprd(Work(I000),Work(I010),T2amp(Ibegin+Ioff),
     +                       Nrow*Ncol)
            Ioff = Ioff + Nrow*Ncol
          Enddo
      Enddo 

      Ibegin = T2ln_aa + T2ln_bb
      Ioff   = 1
      Listw  = 16
      ListD  = 50 

      Do irrep = 1, Nirrep
        Nrow = Irpdpd(Irrep,13)
        Ncol = Irpdpd(Irrep,14)
        I000 = 1
        I010 = I000 + Nrow*Ncol
        Iend = I010 + Nrow*Ncol

        If (Iend .Ge. Maxcor) Call Insmem("get_mbpt_nos",
     +                                     Iend,Maxcor)

        Call Getlst(Work(I000),1,Ncol,2,Irrep,Listw)
        Call Getlst(Work(I010),1,Ncol,2,Irrep,ListD)
        Call vecprd(Work(I000),Work(I010),T2amp(Ibegin+Ioff),
     +              Nrow*Ncol)
        Ioff = Ioff + Nrow*Ncol
      Enddo

#ifdef _DEBUG_LVL0
      Call checksum("T2AA :",T2amp(1),t2ln_aa)
      If (Iuhf .Ne.0) Call checksum("T2BB :",T2amp(1+T2ln_aa),t2ln_bb)
      Call checksum("T2AB :",T2amp(1+T2ln_aa+T2ln_bb),t2ln_ab)
#endif 
   
      Return
      End 
      
      
