      Subroutine Pccd_asymmovvo(Vo,Ov,Tmp1,Tmp2,Lenvo,Nocc,Nvrt)

      Implicit Double Precision(A-H,O-Z)
      Logical Symmetry

      Dimension Vo(Lenvo)
      Dimension Ov(Lenvo)
      Dimension Tmp1(Lenvo)
      Dimension Tmp2(Lenvo)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "sym.com"

      Common /Symm/Symmetry

      Data One,Onem,Ione /1.0D0,-1.0D0,1/

      Irrepx = Ione
      Ispin  = Ione
      Ioff   = Ione

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The VO-MO gradient matrices:@-Entry"
      call output(Vo,1,nvrt,1,nocc,nvrt,nocc,1)
      Write(6,"(a)") "The OV-MO gradient matrices"
      call output(Ov,1,nocc,1,nvrt,nocc,nvrt,1)
#endif
#ifdef _NOSKIP 
      Do Irrepr = 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irrepx)
         Nrow = Vrt(Irrepl,Ispin)
         Ncol = Pop(Irrepr,Ispin)
         Call Dcopy(Ncol*Nrow,Vo(Ioff),1,Tmp1(Ioff),1)
         Call Daxpy(Ncol*Nrow,Onem,Ov(Ioff),1,Vo(Ioff),1)
         Call Daxpy(Ncol*Nrow,Onem,Tmp1(Ioff),1,Ov(Ioff),1)
         Ioff = Ioff + Nrow*Ncol
      Enddo
#endif 

      Nrow = Nocc
      Ncol = Nvrt
      Call Dcopy(Ncol*Nrow,Vo(Ioff),1,Tmp2(Ioff),1)

      Call Transp(Ov,Tmp1,Nvrt,Nocc)
      Call Daxpy(Ncol*Nrow,Onem,Tmp1(Ioff),1,Vo(Ioff),1)

      Call Transp(Tmp2,Tmp1,Nocc,Nvrt)
      Call Daxpy(Ncol*Nrow,Onem,Tmp1(Ioff),1,Ov(Ioff),1)
      
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The VO-MO antisymm. gradient matrices",
     +               "@-Exit"
      call output(Vo,1,nvrt,1,nocc,nvrt,nocc,1)
      Write(6,"(a)") "The OV-MO antisymm. gradient matrices"
      call output(Ov,1,nocc,1,nvrt,nocc,nvrt,1)
      
#endif 

      Return
      End 


