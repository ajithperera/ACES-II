      Subroutine Pccd_antisymmetrize(Htau_pq,Work,Maxcor,Nbas,Nocc,
     +                               Nvrt)

       Implicit Double Precision(A-H,O-Z)
       Logical Sym_packed 

       Dimension Htau_pq(Nbas,Nbas)
       Dimension Work(Maxcor)

#include "machsp.com"
#include "sym.com"
#include "sympop.com"
#include "flags.com"
#include "flags2.com"

      Data Ione, Izero /1,0/

      Irrepx  = Ione
      Length  = Nbas*Nbas

      Sym_packed = .False. 
      If (Sym_packed) Then
         Ndim_oo  = Irpdpd(Irrepx,21)
         Ndim_vv  = Irpdpd(Irrepx,19)
         Ndim_vo  = Irpdpd(Irrepx,9)
      Else
         Ndim_oo = Nocc*Nocc
         Ndim_vv = Nvrt*Nvrt
         Ndim_vo = Nocc*Nvrt
      Endif 

      Lenoo    = Ndim_oo
      Lenvv    = Ndim_vv
      Lenvo    = Ndim_vo

      I000 = Ione
      I010 = I000 + Ndim_oo
      I020 = I010 + Ndim_vv
      I030 = I020 + Ndim_vo
      Iend = I030 + Ndim_vo
      Memleft = Maxcor - Iend
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_antisymmetrize",Iend,
     +                                   Maxcor)

      Call Pccd_frmblocks(Htau_pq,Work(I000),Work(I010),Work(I020),
     +                    Work(I030),Work(Iend),Memleft,Lenoo,Lenvv,
     +                    Lenvo,Nocc,Nvrt,Nbas)

C Store the OV and VO gradients befor antsymmetrize them.      

      Call Putrec(20,"JOBARC","VO_GRADS",Ndim_vo*Iintfp,Work(I020))
      Call Putrec(20,"JOBARC","OV_GRADS",Ndim_vo*Iintfp,Work(I030))

      I040 = Iend
      I050 = I040 + Ndim_oo
      I060 = I050 + Ndim_vv
      I070 = I060 + Ndim_vo 
      Iend = I070 + Ndim_vo 
      Memleft = Maxcor - Iend
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_antisymmetrize",Iend,
     +                                   Maxcor)

      Call Pccd_asymmoovv(Work(I000),Work(I010),Work(I040),Work(I050),
     +                    Lenoo,Lenvv,Nocc,Nvrt)
      Call Pccd_asymmovvo(Work(I020),Work(I030),Work(I060),Work(I070),
     +                    Lenvo,Nocc,Nvrt)

      Call Pccd_putblocks(Htau_pq,Work(I000),Work(I010),Work(I020),
     +                    Work(I030),Work(Iend),Memleft,
     +                    Lenoo,Lenvv,Lenvo,Nocc,Nvrt,Nbas,"Vo_like",
     +                    Sym_packed)

#ifdef _DEBUG_LVLM
      Write(6,"(a)") " The orbital rotation gradient"
      Call output(Htau_pq,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif 
      Return 
      End 


       
