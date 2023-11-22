      Subroutine Pccd_form_nmrcl_grad_hess(Work,Maxcor,Iuhf,Igrad_calc,
     +                                     IHess_calc,Scale)
      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Maxcor)

#include "info.com"
#include "machsp.com"

      Nbasis = Nocco(1) + Nvrto(1)
      Isqlen = Nbasis*Nbasis
      Ndis = 2
      Iplus = Ione
      Imins = Iplus + Isqlen
      Iderv = Imins + Isqlen
      Iend  = Iderv + Isqlen

      Do Idis = 1, Ndis
C  Here we need to Call orbital rotation driver, This code is not 
C complete. Needs lot od work to get numerical gradient and Hessians. 
         
         If (Idis .Eq. Itwo) Scale = -Scale
         Call Pccd_form_htau_pq(Work,Maxcor,Iuhf,Igrad_calc,
     +                          IHess_calc,Scale)
         Call Getrec(20,"JOBARC","OBRTGRDA",Isqlen*IIntfp,Work(Iplus))
         If (Idis .Eq. Itwo) Call Getrec(20,"JOBARC","OBRTGRDA",
     +                                   Isqlen*IIntfp,Work(Imins))

         Dinv = Half*(One/Scale)
         Do I = 1, Isqlen
            Work(Iderv-1+I) = (Work(Iplus-1+I) - Work(Imins-1+I))*Dinv
         Enddo
      Enddo 

      Return 
      End
