      Subroutine Pccd_asymmoovv(Oo,Vv,Oot,Vvt,Lenoo,Lenvv,Nocc,
     +                          Nvrt)

      Implicit Double Precision(A-H,O-Z)

      Dimension Oo(Lenoo)
      Dimension Oot(Lenoo)
      Dimension Vv(Lenvv)
      Dimension Vvt(Lenvv)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "sym.com"
     
      Data Onem,Ione /-1.0D0,1/

      Irrepx = Ione
      Ispin  = Ione
      Ioff   = Ione
      Call Dzero(Vvt,Lenvv)
      Call Dzero(OOt,Lenoo)

      Nrow = Nocc
      Ncol = Nocc
      Call Transp(Oo(Ioff),Oot(Ioff),Ncol,Nrow)
      Call Daxpy(Ncol*Nrow,Onem,Oot(Ioff),1,Oo(Ioff),1)

      Nrow = Nvrt
      Ncol = Nvrt
      Call Transp(Vv(Ioff),Vvt(Ioff),Ncol,Nrow)
      Call Daxpy(Ncol*Nrow,Onem,Vvt(Ioff),1,Vv(Ioff),1)

#ifdef _NOSKIP
      Do Irrepr = 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irrepx)
         Nrow = Pop(Irrepl,Ispin)
         Ncol = Pop(Irrepr,Ispin)
         Call Transp(Oo(Ioff),Oot(Ioff),Ncol,Nrow)
         Call Daxpy(Ncol*Nrow,Onem,Oot(Ioff),1,Oo(Ioff),1)
         Ioff = Ioff + Nrow*Ncol
      Enddo

      Ioff   = Ione
      Do Irrepr = 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irrepx)
	 Nrow = Vrt(Irrepl,Ispin)
         Ncol = Vrt(Irrepr,Ispin)
         Call Transp(Vv(Ioff),Vvt(Ioff),Ncol,Nrow)
         Call Daxpy(Ncol*Nrow,Onem,Vvt(Ioff),1,Vv(Ioff),1)
         Ioff = Ioff + Nrow*Ncol
      Enddo
#endif

#ifdef _DEBUG_LVL0
      write(6,*)
      Write(6,"(a)") "The OO-MO antisymm. gradient matrices"
      call output(OO,1,nocc,1,nocc,nocc,nocc,1)
      Write(6,"(a)") "The VV-MO antisymm. gradient matrices"
      call output(VV,1,nvrt,1,nvrt,nvrt,nvrt,1)
#endif
      Return
      End 


