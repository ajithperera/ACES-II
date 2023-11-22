      Subroutine Pccd_Frmful(Htau,Hoo,Hvv,Hvo,Hov,Work,Maxcor,Nocc,
     +                       Nvrt,Nbas,String,Sym_packed)

       Implicit Double Precision(A-H,O-Z)
       Character*7 String
       Logical Sym_packed 
       Logical Symmetry

       Dimension Htau(Nbas,Nbas)
       Dimension Hoo(Nocc*Nocc)
       Dimension Hvv(Nvrt*Nvrt)
       Dimension Hvo(Nvrt*Nocc)
       Dimension Hov(Nvrt*Nocc)
       Dimension Work(Maxcor)

#include "machsp.com"
#include "sym.com"
#include "sympop.com"
      Common/Symm/Symmetry

      Data Ione, Izero /1,0/

      Irrepx   = Ione
      If (Sym_packed) Then 
          Lenoo    = Irpdpd(Irrepx,21)
          Lenvv    = Irpdpd(Irrepx,19)
          Lenvo    = Irpdpd(Irrepx,9)
      Else
          Lenoo    = Nocc*Nocc
          Lenvv    = Nvrt*Nvrt
          Lenvo    = Nvrt*Nocc
      Endif 

      Call Pccd_putblocks(Htau,Hoo,Hvv,Hvo,Hov,Work,Maxcor,Lenoo,
     +                    Lenvv,Lenvo,Nocc,Nvrt,Nbas,String,
     +                    Sym_packed)

      Return 
      End 


       
