      Subroutine Pccd_hess_update(Hess_old,Grd_old,Hess_vo,Hess_oo,
     +                            Hess_vv,Grd_oo,Grd_vv,Grd_vo,
     +                            Grd_ov,Work,Maxcor,Lhess,Lgrad,
     +                            Nocc,Nvrt,Lenoo,Lenvv,Lenvo,
     +                            Nsize_oo,Nsize_vv,Nsize_vo,
     +                            Icycle)

       Implicit Double Precision (A-H,O-Z)

       Dimension Grd_oo(Lenoo),Grd_vv(Lenvv),Grd_vo(Lenvo)
       Dimension Grd_ov(Lenvo)
       Dimension Hess_old(Lhess),Grd_old(Lgrad)
       Dimension Hess_oo(Nsize_oo),Hess_vv(Nsize_vv),Hess_vo(Nsize_vo)
       Dimension Work(Maxcor)

       Parameter(Tol = 1.0D-15)

#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"

      Common /Symm/Symmetry

      Data Ione,Izero,One,Half,Dnull,Two,Onem/1,0,1.0D0,0.50D0,
     +                                        0.0D0,2.0D0,-1.0D0/

      Len_grad    = Lenoo
      Call Getrec(20,"JOBARC","OO_ROT_H",Nsize_oo*Iintfp,
     +            Hess_old(Ione))
      Call Getrec(20,"JOBARC","OO_ROT_G",Lenoo*Iintfp,Grd_old(Ione))

      I000 = Ione
      I010 = I000 + Len_grad
      Call Dcopy(Lenoo,Grd_old,1,Work(I000),1)

      Grad_norm_o = Ddot(Lenoo,Grd_old(I000),1,Grd_old(I000),1)
      Grad_norm_n = Ddot(Lenoo,Work(I000),1,Work(I000),1)

      If (Grad_norm_o .Gt. Tol .and. Grad_norm_n .Gt. Tol) Then 
         I020 = I010 + 7*Len_grad*Len_grad
         I030 = I020 + 3*Len_grad*Len_grad
         Iend = I030 + Len_grad
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_bfgs",Iend,Maxcor)

         Call Getrec(20,"JOBARC","OO_HG_ST",Lenoo*Iintfp,Work(I030))
         Call Pccd_bfgs(Grd_oo,Hess_old,Work(I000),Work(I030),
     +                  Work(I020),Len_grad)
         Call Dcopy(Nsize_oo,Hess_old,1,Hess_oo,1)
      Else
         Call Dcopy(Nsize_oo,Hess_old,1,Hess_oo,1)
      Endif 

#ifdef _DEBUG_LVL0
      write(6,"(a)") " @-pccd_hess_update: updated OO Hessian"
      call output(Hess_oo,1,Lenoo,1,Lenoo,Lenoo,Lenoo,1)
#endif 

      Len_grad    = Lenvo
      Call Getrec(20,"JOBARC","VO_ROT_H",Nsize_vo*Iintfp,
     +            Hess_old(Ione))
      Call Getrec(20,"JOBARC","VO_ROT_G",Lenvo*Iintfp,
     +            Grd_old(Ione))
      I000 = Ione
      I010 = I000 + Len_grad
      Call Dcopy(Lenvo,Grd_old,1,Work(I000),1)

      Grad_norm_o = Ddot(Lenvo,Grd_old(I000),1,Grd_old(I000),1)
      Grad_norm_n = Ddot(Lenvo,Work(I000),1,Work(I000),1)

      If (Grad_norm_o .Gt. Tol .and. Grad_norm_n .Gt. Tol) Then 
         I020 = I010 + 7*Len_grad*Len_grad
         I030 = I020 + 3*Len_grad*Len_grad
         Iend = I030 + Len_grad
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_bfgs",Iend,Maxcor)

         Call Getrec(20,"JOBARC","VO_HG_ST",Lenvo*Iintfp,Work(I030))
         Call Pccd_bfgs(Grd_vo,Hess_old,Work(I000),Work(I030),
     +                  Work(I020),Len_grad)
         Call Dcopy(Nsize_vo,Hess_old,1,Hess_vo,1)
      Else
         Call Dcopy(Nsize_vo,Hess_old,1,Hess_vo,1)
      Endif 

#ifdef _DEBUG_LVL0
      write(6,"(a)") "@-pccd_hess_update: Updated VO Hessian"
      call output(Hess_vo,1,Lenvo,1,Lenvo,Lenvo,Lenvo,1)
#endif 

      Len_grad    = Lenvv
      Call Getrec(20,"JOBARC","VV_ROT_H",Nsize_vv*Iintfp,
     +            Hess_old(Ione))
      Call Getrec(20,"JOBARC","VV_ROT_G",Lenvv*Iintfp,
     +            Grd_old(Ione))
      I000 = Ione
      I010 = I000 + Len_grad
      Call Dcopy(Lenvv,Grd_old,1,Work(I000),1)

      Grad_norm_o = Ddot(Lenvv,Grd_old(I000),1,Grd_old(I000),1)
      Grad_norm_n = Ddot(Lenvv,Work(I000),1,Work(I000),1)

      If (Grad_norm_o .Gt. Tol .and. Grad_norm_n .Gt. Tol) Then 
         I020 = I010 + 7*Len_grad*Len_grad
         I030 = I020 + 3*Len_grad*Len_grad
         Iend = I030 + Len_grad
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_bfgs",Iend,Maxcor)

         Call Getrec(20,"JOBARC","VV_HG_ST",Lenvv*Iintfp,Work(I030))
         Call Pccd_bfgs(Grd_vv,Hess_old,Work(I000),Work(I030),
     +                  Work(I020),Len_grad)
         Call Dcopy(Nsize_vv,Hess_old,1,Hess_vv,1)
      Else
         Call Dcopy(Nsize_vv,Hess_old,1,Hess_vv,1)
      Endif 
#ifdef _DEBUG_LVL0
      write(6,"(a)") "@-pccd_hess_update: Updated VV Hessian"
      call output(Hess_vv,1,Lenvv,1,Lenvv,Lenvv,Lenvv,1)
#endif 

      Return
      End
   
 
