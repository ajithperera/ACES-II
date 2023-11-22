#include "flags.h"
      Subroutine Tdcc_form_mu_0_s(Work,Memleft,Doo,Dvv,Dvo,Mu_0_S,
     +                            Irrepx,Lenoo,Lenvv,Lenvo,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)
      Dimension Doo(Lenoo),Dvv(Lenvv),Dvo(Lenvo)

      Double Precision Mu_0_s(Lenvo)
      Double Precision Mone
      Integer dissiz_aa,dissiz_bb,dissiz_ab
 
      Data Zero, One, Mone, Two /0.0D0, 1.0D0, -1.0D0, 2.0D0/

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "symloc.com"
#include "sym.com"
#ifdef _DEBUG_LVLM
        Write(6,*)
        Write(6,"(a)") "---Entered tdcc_form_mu_0_S---"
        Write(6,*)
#endif
      Mu0_ph_list = 390

      Do Ispin = 1, 1+Iuhf 
         
         Ivo_0 = 1 + (Ispin-1) * Irpdpd(Irrepx,9)

         Do Irrep_i = 1, Nirrep
            Irrep_a = Dirprd(Irrep_i,Irrepx)

            Ndim_i = Pop(Irrep_i,Ispin)
            Ndim_a = Vrt(Irrep_a,Ispin)
               
            Ivo_z = Ivo_0 + (Isymoff(Irrep_i,Irrepx,8+Ispin)-1)
            Ivo_d = Ivo_z

C z(A,I) = Mu(A,I) (and Z(a,i) = Mu(a,i))

            Ndim_ai = Ndim_i * Ndim_a
            Call Dcopy(Ndim_ai,Dvo(Ivo_d),1,Mu_0_S(Ivo_z),1)

         Enddo

         Call Putlst(Mu_0_S(Ivo_0),1,1,1,Ispin,Mu0_ph_list)

      Enddo

#ifdef _DEBUG_LVL0
      do ispin=1, 1+iuhf
      Call Getlst(Mu_0_S,1,1,1,Ispin,Mu0_ph_list)
      ioff = 1 + (Ispin-1)*Irpdpd(Irrepx,9)
      do irrep_i=1,Nirrep
         irrep_a=dirprd(irrep_i,Irrepx)
         ndim_i = Pop(irrep_i,ispin)
         Ndim_a = Vrt(irrep_a,ispin)
         Write(6,*)
         write(6,"(a,a,5(1x,i2))") "irrep_i,Irrep_a,Ndim_i,",
     +                     "Ndim_a,Ioff: ",irrep_i,Irrep_a,
     +                      Ndim_i,Ndim_a,Ioff
         call output(Mu_0_S(ioff),1,ndim_a,1,Ndim_i,Ndim_a,
     +                Ndim_i,1)
         Ioff = Ioff + ndim_i*Ndim_a
      enddo
      enddo
#endif

      Return
      End
     
           
        
