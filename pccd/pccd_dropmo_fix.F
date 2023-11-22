      Subroutine  Pccd_dropmo_fix(Grad,Grad_t,Nbas,Ndrop)

      Implicit Double Precision(A-H,O-Z)
      Dimension Grad(Nbas,NBas), Idrop(Ndrop), Grad_t(Nbas,Nbas)

#include "syminf.com"
#include "sym.com"
#include "machsp.com"
#include "sympop.com"
#include "info.com"

      Data Azero /0.0D0/

      Call Getrec(20,"JOBARC","MODROPA",Ndrop,Idrop)

      Noc_drp = 0
      Nvr_drp = 0

      Do I = 1, Ndrop
         If (Idrop(i) .Le. Nocco(1)) Then
             Noc_drp = Noc_drp + 1
         Else
             Nvr_drp = Nvr_drp +1 
         Endif
      Enddo 

C We always assume that the first Noc_drp orbitals and last Nbas-Nvr_drp
C orbitals are dropped. Zero out the corresponding blocks of the
C gradient matrix

       Do I = 1, Noc_drp
          Do J = 1, Nbas 
             Grad(J,I) = Azero
          Enddo
       Enddo  

       Do I = 1, Nbas
          Do J = 1, Noc_drp 
             Grad(J,I) = Azero
          Enddo
       Enddo 


       Nstart = Nbas-Nvr_drp+1
       Do I = Nstart, Nbas
          Do J = 1, Nbas 
             Grad(J,I) = Azero
          Enddo
       Enddo  

       Do I = 1, Nbas
          Do J = Nstart, Nbas 
             Grad(J,I) = Azero
          Enddo
       Enddo  

       Call Dzero(Grad_t,Nbas*Nbas)

       K = 0
       L = 0
       Nstart = Nbas-Nvr_drp
       Do I = Noc_drp+1, Nstart
             K = K + 1
             L = 0
          Do J = Noc_drp+1, Nstart
             L = L + 1
             Grad_t(L,K) =Grad(J,I)
             Print*, "J and I", J,I,Grad(J,I)
             Print*, "K and L", K,L,Grad_t(L,K)
          Enddo
       ENddo
       Print*,"Test1"
       Call output(Grad,1,Nbas,1,Nbas,Nbas,Nbas,1)
   
       Ndrop_mo = Nbas - Ndrop 
       Call output(Grad_t,1,Ndrop_mo,1,Ndrop_mo,nbas,nbas1)
       Call Dcopy(Ndrop_mo*Ndrop_mo,Grad_t,1,Grad,1)
       Print*,"Test2"
       Call output(Grad,1,Ndrop_mo,1,Ndrop_mo,Nbas,Nbas,1)
      
       Return 
       End 
