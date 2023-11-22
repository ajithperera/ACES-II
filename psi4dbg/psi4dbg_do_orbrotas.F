#include "flags.h"
      Subroutine Psi4dbg_do_orbrots(Work,Maxcor,Iuhf,IGrad_calc,
     +                              IHess_calc,Scale,OOmicroItr,
     +                              OOmacroItr,OOtotalItr)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxcor)

#include "info.com"
#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "flags2.com"

      Data Ione,Itwo,One,Half /1,2,1.0D0,0.50D0/

      Call Psi4dbg_form_htau_pq(Work,Maxcor,Iuhf,Igrad_calc,
     +                          IHess_calc,Scale)

      Call GetSDInfo(nocc,nvirt,nbas,nvecDim,Naobfns)


      If (Iflags2(h_IFLAGS2_oo_method) .EQ. 0) then ! BB Steep. Desc.
         call Drive_BB_SD(Work,Maxcor,Iuhf,OOmicroItr,OOmacroItr,
     +                    nocc,nvirt,nvecDim)
      else if (Iflags2(h_IFLAGS2_oo_method) .EQ. 1) then !BFGS
        print*,'Have not include L-BFGS yet; to do'
      else if (Iflags2(h_IFLAGS2_oo_method) .EQ. 2) then !FullNR
        print*,'Have not included FullNR yet; to do'
      else if (Iflags2(h_IFLAGS2_oo_method) .EQ. 3) then !AMSgrad
          call Drive_AMSgrad(Work,Maxcor,Iuhf,OOmicroItr,OOmacroItr,
     +                        nocc,nvirt,nvecDim)
      else if (Iflags2(h_IFLAGS2_oo_method) .EQ. 4) then !AMSgrad
          Call Psi4dbg_rotgrd(Work,Maxcor,Nbas,Nocc,Nvirt)
      endif

      IFLAGS(h_IFLAGS_scf_maxcyc)=0 
C      IFLAGS(h_IFLAGS_non_hf)=1 
      CALL PUTREC(20,'JOBARC','IFLAGS  ',100,IFLAGS)
      Return
      End 
 
