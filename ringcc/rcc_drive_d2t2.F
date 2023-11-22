      Subroutine Rcc_drive_d2t2(Work,Maxcor,Iuhf)

      Implicit Double Precision (A-H, O-Z)
      Dimension Work(Maxcor)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "flags.com"
#include "info.com"

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "--------Entering rcc_drive_d2t2---------"
      Write(6,*)
#endif
      
C Obtain the maximum length of the Denominator AA and BB array

      Maxaabb = 1
      Do Ispin = 1, Iuhf+1
         Do Irrep = 1, Nirrep
            Nrow = Irpdpd(Irrep,4+Ispin)
            NCol = Irpdpd(Irrep,6+Ispin)
            Length = Nrow * Ncol
            If (Length .LT. 0) Call Trap_intovf("Rcc_make_d2t2",1)
            Maxaabb=Max(Maxaabb,Length)
         Enddo 
      Enddo 

      Nbasis = Nocco(1) + Nvrto(1)
      Itmp1  = Nocco(1) * Nocco(1)
      Itmp1  = Nvrto(2) * Nvrto(2)
      Isize  = Itmp1 + Itmp2 + 2*Nirrep 
C  
      I000 = 1
      I010 = I000 + 2*Nbasis
      I020 = I010 + Maxaabb
      I030 = I020 + Isize
      I030 = Iend 
      If (Iend .LT. 0) Call Trap_intovf("Rcc_make_d2t2",2)
      If (Iend .GT. Maxcor) Call Insmem("Rcc_make_d2t2",Iend,Maxcor)
      Call Rcc_make_d2(Work(I000),Work(I010),Work(I020),Nbasis,
     +                 Isize,Maxaabb,Iuhf)

      Call Rcc_make_t2(Work,Maxcor,Iuhf)
     
      Return
      End 
