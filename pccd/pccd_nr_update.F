      Subroutine Pccd_nr_update(Hess_vo,Hess_oo,Hess_vv,Grd_oo,Grd_vv,
     +                          Grd_vo,Grd_ov,Work,Maxcor,Lenoo,
     +                          Lenvv,Lenvo,Nocc,Nvrt,Wrt2jarc)

       Implicit Double Precision (A-H,O-Z)
       Logical Wrt2jarc 

       Dimension Grd_oo(Lenoo),Grd_vv(Lenvv),Grd_vo(Lenvo)
       Dimension Grd_ov(Lenvo)
       Dimension Hess_oo(Lenoo,Lenoo),Hess_vo(Lenvo,Lenvo)
       Dimension Hess_vv(Lenvv,Lenvv)
       Dimension Work(Maxcor)
       
#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"

      Common /Symm/Symmetry

      Data Ione,Izero,One,Half,Dnull,Two,Onem/1,0,1.0D0,0.50D0,
     +                                        0.0D0,2.0D0,-1.0D0/
     
      I000 = Ione
      I010 = I000 + Max(Lenvv,Lenoo,Lenvo)
      Iend = I010 + Lenvo*Lenvo

      Call Dgemm("N","N",Lenvo,Lenvo,Lenvo,Onem,Hess_vo,Lenvo,Grd_ov,
     +            Lenvo,Dnull,Work,Lenvo)
      Call Dcopy(Lenvo,Work,1,Grd_ov,1)
    
      Call Dgemm("T","N",Lenvo,Lenvo,Lenvo,Onem,Hess_vo,Lenvo,Grd_vo,
     +            Lenvo,Dnull,Work,Lenvo)
      Call Dcopy(Lenvo,Work,1,Grd_vo,1)

      Call Dgemm("N","N",Lenoo,Lenoo,Lenoo,Onem,Hess_oo,Lenoo,Grd_oo,
     +            Lenoo, Dnull,Work,Lenoo)
      Call Dcopy(Lenoo,Work,1,Grd_oo,1)


      Call Dgemm("N","N",Lenvv,Lenvv,Lenvv,Onem,Hess_vv,Lenvv,Grd_vv,
     +            Lenvv,Dnull,Work,Lenvv)
      Call Dcopy(Lenvv,Work,1,Grd_vv,1)

#ifdef _DEBUG_LVL0
      Write(6,*) "-----------Printing begins @-nr_update------"
       Write(6,*)
       Write(6,"(a)") "The OO:H(^-1)G gradients"
       call output(Grd_oo,1,Nocc,1,Nocc,Nocc,Nocc,1)
       Write(6,"(a)") "The VV:H(^-1)G gradients"
       call output(Grd_vv,1,Nvrt,1,Nvrt,Nvrt,Nocc,1)
       Write(6,"(a)") "The VO:H(^-1)G gradients"
       call output(Grd_vo,1,Nvrt,1,Nocc,Nvrt,Nocc,1)
       Write(6,"(a)") "The OV:H(^-1)G gradients"
       call output(Grd_ov,1,Nocc,1,Nvrt,Nocc,Nvrt,1)
      Write(6,*) "-----------Printing ends @-nr_update--------"
#endif
      If (Wrt2jarc) Then
         Call Putrec(20,"JOBARC","OO_HG_ST",Lenoo*Iintfp,Grd_oo)
         Call Putrec(20,"JOBARC","OV_HG_ST",Lenvo*Iintfp,Grd_vo)
         Call Putrec(20,"JOBARC","VO_HG_ST",Lenvo*Iintfp,Grd_ov)
         Call Putrec(20,"JOBARC","VV_HG_ST",Lenvv*Iintfp,Grd_vv)
      Endif 

      Return
      End
   
 
